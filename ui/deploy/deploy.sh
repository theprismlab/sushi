#!/usr/bin/env bash
# Deploy one environment of the sushi UI. Idempotent: exits immediately if the
# branch has not moved, so it is safe to run from cron every couple of minutes.
#
#   deploy.sh develop            # deploy if origin/develop moved
#   deploy.sh main --force       # redeploy even if nothing moved
#
# Run as root: it restarts the systemd unit. The service itself runs as jenkins.

# The outer brace matters: bash reads to the closing brace before executing
# anything, so a deploy that updates this very file cannot corrupt the run in
# progress.
{
set -euo pipefail

BRANCH="${1:?usage: deploy.sh <develop|main> [--force]}"
FORCE="${2:-}"

ROOT="${SUSHI_UI_ROOT:-/opt/sushi-ui}"
CHECKOUT="$ROOT/$BRANCH"
UNIT="sushi-ui@$BRANCH"
NODE_IMAGE="${NODE_IMAGE:-docker.io/library/node:20-alpine}"

# cron gives a minimal PATH; uv, podman and systemctl all need to be findable.
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:${PATH:-}"

log() { printf '%s [%s] %s\n' "$(date -u +%FT%TZ)" "$BRANCH" "$*"; }
die() { log "ERROR: $*"; exit 1; }

case "$BRANCH" in
  develop|main) ;;
  *) die "refusing to deploy branch '$BRANCH'; only develop and main have environments" ;;
esac

[[ -d "$CHECKOUT/.git" ]] || die "no git checkout at $CHECKOUT -- run install.sh first"

# One deploy per environment at a time. npm ci is slow enough that overlapping
# cron ticks are a real possibility.
exec 9>"/var/lock/sushi-ui-deploy.$BRANCH"
flock -n 9 || { log "another deploy is already running; skipping"; exit 0; }

cd "$CHECKOUT"

PORT="$(sed -n 's/^PORT=//p' "$CHECKOUT/ui/deploy/env.$BRANCH" | head -1)"
[[ -n "$PORT" ]] || die "no PORT in ui/deploy/env.$BRANCH"

git fetch --quiet origin "$BRANCH"
CURRENT="$(git rev-parse HEAD)"
TARGET="$(git rev-parse "origin/$BRANCH")"

if [[ "$CURRENT" == "$TARGET" && "$FORCE" != "--force" ]]; then
  exit 0   # nothing to do; stay quiet so the cron log stays readable
fi

build_and_start() {
  local commit="$1"
  git reset --hard --quiet "$commit"

  # The VM has no node, so the SPA is built in a container. node_modules and
  # dist are gitignored, so reset --hard leaves the npm cache in place.
  podman run --rm \
    -v "$CHECKOUT/ui/frontend:/app:Z" \
    -w /app "$NODE_IMAGE" \
    sh -c 'npm ci --no-audit --no-fund && npm run build'

  uv sync --directory "$CHECKOUT/ui/backend" --frozen --no-dev

  systemctl restart "$UNIT"
}

healthy() {
  local i
  for i in $(seq 1 30); do
    if curl -fsS --max-time 3 "http://localhost:$PORT/api/health" >/dev/null 2>&1; then
      return 0
    fi
    sleep 1
  done
  return 1
}

log "deploying ${CURRENT:0:8} -> ${TARGET:0:8}"
build_and_start "$TARGET"

if healthy; then
  log "healthy on :$PORT at ${TARGET:0:8}"
  exit 0
fi

# An unattended deploy must not be able to leave an environment dead.
log "health check failed after 30s; rolling back to ${CURRENT:0:8}"
if build_and_start "$CURRENT" && healthy; then
  log "rolled back to ${CURRENT:0:8}; ${TARGET:0:8} was NOT deployed"
  exit 1
fi

die "rollback to ${CURRENT:0:8} also failed -- $UNIT is down, needs a human"
}
exit
