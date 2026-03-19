#!/usr/bin/env bash
set -e

echo "PWD=$(pwd)"
echo "Contents of repo root:"
ls -1

# Build, using the repo root (.) as the context,
# and pointing Docker at the Dockerfile under ./docker/
podman build \
  --pull \
  --no-cache \
  --build-arg GIT_COMMIT="$(git rev-parse HEAD)" \
  -f docker/Dockerfile \
  -t sushi-podman:production \
  .
