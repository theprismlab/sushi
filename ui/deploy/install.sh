#!/usr/bin/env bash
# One-time bootstrap on vercingetorix-r8. Run as root. Idempotent -- re-run it
# after changing the unit file, after adding an environment, or once `main` has
# caught up with the UI code.
#
#   sudo ./install.sh                  # provision every environment that is ready
#   sudo ./install.sh develop          # just this one

set -euo pipefail

REPO="${REPO:-git@github.com:theprismlab/sushi.git}"
ROOT="${SUSHI_UI_ROOT:-/opt/sushi-ui}"
DB_DIR="${DB_DIR:-/local/jenkins}"
SERVICE_USER="${SERVICE_USER:-jenkins}"

# Read the unit template from this script's own directory rather than from the
# main checkout: main does not necessarily have the UI code yet.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

BRANCHES=("$@")
[[ ${#BRANCHES[@]} -eq 0 ]] && BRANCHES=(main develop)

[[ $EUID -eq 0 ]] || { echo "run as root"; exit 1; }
command -v uv >/dev/null || { echo "uv is not on PATH"; exit 1; }
command -v podman >/dev/null || { echo "podman is not on PATH"; exit 1; }

install -m644 "$SCRIPT_DIR/sushi-ui@.service" /etc/systemd/system/
systemctl daemon-reload

# The run databases are created by the app on first start; the directory just
# has to be writable by the service user.
install -d -o "$SERVICE_USER" -g "$SERVICE_USER" "$DB_DIR"

ready=()
for branch in "${BRANCHES[@]}"; do
  if [[ ! -d "$ROOT/$branch/.git" ]]; then
    echo "cloning $branch into $ROOT/$branch"
    git clone --quiet --branch "$branch" "$REPO" "$ROOT/$branch"
  else
    git -C "$ROOT/$branch" fetch --quiet origin "$branch"
    git -C "$ROOT/$branch" reset --hard --quiet "origin/$branch"
  fi

  # A branch that has not merged the UI yet has no environment to provision.
  # Skip it rather than writing a cron entry that fails every two minutes.
  if [[ -f "$ROOT/$branch/ui/deploy/env.$branch" ]]; then
    ready+=("$branch")
  else
    echo "SKIP $branch: no ui/deploy/env.$branch on that branch yet."
    echo "     Merge the UI into $branch, then re-run this script."
  fi
done

[[ ${#ready[@]} -gt 0 ]] || { echo "nothing to provision"; exit 1; }

# Poll for pushes. Deliberately not a GitHub webhook or a hosted CI runner:
# this VM is on an internal network that github.com cannot reach.
{
  echo "# Auto-deploy the sushi UI. Managed by ui/deploy/install.sh -- edits are"
  echo "# overwritten on the next run. deploy.sh exits silently when the branch"
  echo "# has not moved, so this is cheap."
  echo "PATH=/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"
  for branch in "${ready[@]}"; do
    echo "*/2 * * * * root $ROOT/$branch/ui/deploy/deploy.sh $branch >> /var/log/sushi-ui-deploy.log 2>&1"
  done
} > /etc/cron.d/sushi-ui
chmod 644 /etc/cron.d/sushi-ui

cat > /etc/logrotate.d/sushi-ui <<'ROTATE'
/var/log/sushi-ui-deploy.log {
    weekly
    rotate 8
    compress
    missingok
    notifempty
}
ROTATE

for branch in "${ready[@]}"; do
  "$ROOT/$branch/ui/deploy/deploy.sh" "$branch" --force
  systemctl enable "sushi-ui@$branch"
  port="$(sed -n 's/^PORT=//p' "$ROOT/$branch/ui/deploy/env.$branch" | head -1)"
  echo "$branch -> http://$(hostname):$port"
done

echo
echo "These ports must stay on the internal network -- the UI has no auth."
