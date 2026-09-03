# Build context is ui/, not the repo root.
#
#   podman build -t sushi-ui:develop -f ui/Dockerfile ui
#
# One image holds the built SPA and the API, so a deploy is a single
# stop/run and there is nothing to install on the host: no node, no venv.

# ---- Stage 1: the SPA. The only place node exists; the VM has none. --------
FROM docker.io/library/node:20-alpine AS frontend
WORKDIR /app
# Copy the manifests alone first so npm ci is cached until deps actually change.
COPY frontend/package.json frontend/package-lock.json ./
RUN npm ci --no-audit --no-fund
COPY frontend/ ./
RUN npm run build

# ---- Stage 2: the service. -------------------------------------------------
FROM docker.io/library/python:3.11-slim
COPY --from=ghcr.io/astral-sh/uv:0.5.11 /uv /usr/local/bin/uv

# The image already has the interpreter uv.lock resolves against; never let uv
# go and fetch another one.
ENV UV_PYTHON_DOWNLOADS=never \
    UV_LINK_MODE=copy \
    PYTHONUNBUFFERED=1

WORKDIR /app/backend
COPY backend/pyproject.toml backend/uv.lock ./
RUN uv sync --frozen --no-dev

COPY backend/ ./
# app.py serves the SPA from ../frontend/dist relative to itself.
COPY --from=frontend /app/dist /app/frontend/dist

ENV PATH="/app/backend/.venv/bin:$PATH"

# PORT differs per environment, so it has to be read at runtime.
CMD ["sh", "-c", "exec uvicorn app:app --host 0.0.0.0 --port ${PORT:-8100}"]
