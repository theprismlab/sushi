#!/usr/bin/env bash
# One-time bootstrap on vercingetorix-r8. Run as root. Idempotent -- safe to
# re-run after changing the unit file or the cron schedule.
#
#   sudo REPO=git@github.com:cmap/sushi.git ./install.sh

set -euo pipefail

REPO="${REPO:-https://github.com/cmap/sushi.git}"
ROOT="${SUSHI_UI_ROOT:-/opt/sushi-ui}"
DB_DIR="${DB_DIR:-/local/jenkins}"
SERVICE_USER="${SERVICE_USER:-jenkins}"

[[ $EUID -eq 0 ]] || { echo "run as root"; exit 1; }
command -v uv >/dev/null || { echo "uv is not on PATH"; exit 1; }
command -v podman >/dev/null || { echo "podman is not on PATH"; exit 1; }

for branch in main develop; do
  if [[ ! -d "$ROOT/$branch/.git" ]]; then
    echo "cloning $branch into $ROOT/$branch"
    git clone --branch "$branch" "$REPO" "$ROOT/$branch"
  fi
done

install -m644 "$ROOT/main/ui/deploy/sushi-ui@.service" /etc/systemd/system/
systemctl daemon-reload

# The run databases are created by the app on first start; make sure the
# directory is writable by the service user.
install -d -o "$SERVICE_USER" -g "$SERVICE_USER" "$DB_DIR"

# Poll for pushes. Deliberately not a GitHub webhook or a hosted CI runner:
# this VM is on an internal network that github.com cannot reach.
cat > /etc/cron.d/sushi-ui <<CRON
# Auto-deploy the sushi UI. Both scripts exit immediately when the branch
# has not moved, so this is cheap.
PATH=/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
*/2 * * * * root $ROOT/develop/ui/deploy/deploy.sh develop >> /var/log/sushi-ui-deploy.log 2>&1
*/2 * * * * root $ROOT/main/ui/deploy/deploy.sh main       >> /var/log/sushi-ui-deploy.log 2>&1
CRON
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

# First deploy builds the SPA and the venv, then enables the units.
"$ROOT/main/ui/deploy/deploy.sh" main --force
"$ROOT/develop/ui/deploy/deploy.sh" develop --force
systemctl enable sushi-ui@main sushi-ui@develop

echo
echo "production : http://$(hostname):8100"
echo "develop    : http://$(hostname):8101"
echo
echo "Both ports must stay on the internal network -- the UI has no auth."
