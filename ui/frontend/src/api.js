async function request(path, options = {}) {
  const res = await fetch(`/api${path}`, {
    headers: { 'Content-Type': 'application/json' },
    credentials: 'same-origin', // carry the session cookie
    ...options,
  })
  const body = res.headers.get('content-type')?.includes('json') ? await res.json() : await res.text()
  if (!res.ok) {
    const detail = body?.detail ?? body
    throw Object.assign(new Error(typeof detail === 'string' ? detail : 'Request failed'), { detail })
  }
  return body
}

export const api = {
  catalog: () => request('/catalog'),
  health: () => request('/health'),
  session: () => request('/session'),
  signIn: (user, token) =>
    request('/session', { method: 'POST', body: JSON.stringify({ user, token }) }),
  signOut: () => request('/session', { method: 'DELETE' }),
  buildPaths: () => request('/build-paths'),
  ls: (path) => request(`/ls?path=${encodeURIComponent(path || '')}`),
  suggestions: () => request('/suggestions'),
  screenTypes: () => request('/screen-types'),
  saveScreenTypes: (presets) =>
    request('/screen-types', { method: 'PUT', body: JSON.stringify({ presets }) }),
  resetScreenTypes: () => request('/screen-types/reset', { method: 'POST' }),
  pathDefaults: (path, preset) =>
    request(`/path-defaults?path=${encodeURIComponent(path || '')}&preset=${encodeURIComponent(preset || '')}`),
  preflight: (body) => request('/preflight', { method: 'POST', body: JSON.stringify(body) }),
  launch: (body) => request('/runs', { method: 'POST', body: JSON.stringify(body) }),
  runs: (params) => request(`/runs?${new URLSearchParams(params)}`),
  run: (id) => request(`/runs/${id}`),
  log: (id, start) => request(`/runs/${id}/log?start=${start}`),
  outputs: (id) => request(`/runs/${id}/outputs`),
  setNotes: (id, notes) => request(`/runs/${id}`, { method: 'PATCH', body: JSON.stringify({ notes }) }),
  stop: (id) => request(`/runs/${id}/stop`, { method: 'POST' }),
}

export const RUNNING_STATES = ['QUEUED', 'RUNNING']

export function formatBytes(n) {
  if (n == null) return '—'
  const units = ['B', 'KB', 'MB', 'GB', 'TB']
  let i = 0
  while (n >= 1024 && i < units.length - 1) { n /= 1024; i += 1 }
  return `${n < 10 && i > 0 ? n.toFixed(1) : Math.round(n)} ${units[i]}`
}

export function formatDuration(ms) {
  if (!ms) return '—'
  const s = Math.round(ms / 1000)
  if (s < 60) return `${s}s`
  const m = Math.floor(s / 60)
  if (m < 60) return `${m}m ${s % 60}s`
  return `${Math.floor(m / 60)}h ${m % 60}m`
}

export function formatTime(iso) {
  if (!iso) return '—'
  return new Date(iso).toLocaleString(undefined, {
    year: 'numeric', month: 'short', day: 'numeric', hour: '2-digit', minute: '2-digit',
  })
}
