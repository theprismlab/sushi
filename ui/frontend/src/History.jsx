import { useEffect, useState } from 'react'
import { Link, useSearchParams } from 'react-router-dom'
import { api, formatDuration, formatTime, RUNNING_STATES } from './api.js'
import { useCatalog } from './App.jsx'

const PAGE = 50
const STATUSES = ['QUEUED', 'RUNNING', 'SUCCESS', 'FAILURE', 'ABORTED', 'UNSTABLE', 'ERROR']

export default function History() {
  const catalog = useCatalog()
  const [params, setParams] = useSearchParams()
  const [data, setData] = useState(null)
  const [error, setError] = useState(null)

  const filters = {
    preset: params.get('preset') || '',
    status: params.get('status') || '',
    build_name: params.get('build_name') || '',
    launched_by: params.get('launched_by') || '',
  }
  const page = Number(params.get('page') || 0)

  function setFilter(key, value) {
    const next = new URLSearchParams(params)
    value ? next.set(key, value) : next.delete(key)
    next.delete('page') // a changed filter invalidates the current page
    setParams(next)
  }

  function setPage(value) {
    const next = new URLSearchParams(params)
    value > 0 ? next.set('page', String(value)) : next.delete('page')
    setParams(next)
  }

  useEffect(() => {
    let cancelled = false
    const load = () => {
      const query = { limit: PAGE, offset: page * PAGE }
      Object.entries(filters).forEach(([k, v]) => { if (v) query[k] = v })
      api.runs(query)
        .then((d) => { if (!cancelled) { setData(d); setError(null) } })
        .catch((e) => { if (!cancelled) setError(e.message) })
    }
    load()
    // Reading the list is what refreshes statuses from Jenkins, so poll while
    // anything is in flight.
    const timer = setInterval(load, 15000)
    return () => { cancelled = true; clearInterval(timer) }
  }, [params.toString()])

  if (error) return <div className="banner error">{error}</div>
  if (!data) return <div className="muted">Loading runs…</div>

  return (
    <div className="history">
      <div className="filters card">
        <label>
          <span className="muted">Screen type</span>
          <select value={filters.preset} onChange={(e) => setFilter('preset', e.target.value)}>
            <option value="">all</option>
            {catalog.presets.map((p) => <option key={p.id} value={p.id}>{p.id}</option>)}
          </select>
        </label>
        <label>
          <span className="muted">Status</span>
          <select value={filters.status} onChange={(e) => setFilter('status', e.target.value)}>
            <option value="">all</option>
            {STATUSES.map((s) => <option key={s} value={s}>{s}</option>)}
          </select>
        </label>
        <label>
          <span className="muted">Build name</span>
          <input value={filters.build_name} onChange={(e) => setFilter('build_name', e.target.value)}
                 placeholder="contains…" />
        </label>
        <label>
          <span className="muted">Launched by</span>
          <input value={filters.launched_by} onChange={(e) => setFilter('launched_by', e.target.value)} />
        </label>
        <span className="muted count">{data.total} runs</span>
      </div>

      {data.runs.length === 0 ? (
        <div className="card muted">
          No runs recorded yet. Builds launched from the Jenkins form directly will not appear here.
        </div>
      ) : (
        <table className="runs card">
          <thead>
            <tr>
              <th>Status</th><th>Build</th><th>Type</th><th>Launched by</th>
              <th>Started</th><th>Duration</th><th>Modules</th><th>Log</th>
            </tr>
          </thead>
          <tbody>
            {data.runs.map((run) => (
              <tr key={run.id}>
                <td><Status value={run.status} /></td>
                <td>
                  <Link to={`/runs/${run.id}`}><strong>{run.build_name}</strong></Link>
                  <div className="muted small">#{run.id}{run.build_number ? ` · jenkins #${run.build_number}` : ''}</div>
                </td>
                <td>
                  {run.preset}
                  {run.screen_type !== run.preset && (
                    <div className="muted small" title="SCREEN_TYPE parameter differs from the screen type directory">
                      as {run.screen_type}
                    </div>
                  )}
                </td>
                <td>{run.launched_by}</td>
                <td className="small">{formatTime(run.started_at || run.created_at)}</td>
                <td className="small">{formatDuration(run.duration_ms)}</td>
                <td className="small">{run.enabled_modules.length}</td>
                <td className="small">
                  {run.jenkins_url && <a href={run.jenkins_url} target="_blank" rel="noreferrer">jenkins</a>}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      )}

      {data.total > PAGE && (
        <div className="pager">
          <button disabled={page === 0} onClick={() => setPage(page - 1)}>← newer</button>
          <span className="muted">{page * PAGE + 1}–{Math.min((page + 1) * PAGE, data.total)} of {data.total}</span>
          <button disabled={(page + 1) * PAGE >= data.total}
                  onClick={() => setPage(page + 1)}>older →</button>
        </div>
      )}
    </div>
  )
}

export function Status({ value }) {
  const kind = value === 'SUCCESS' ? 'ok'
    : RUNNING_STATES.includes(value) ? 'busy'
    : value === 'UNSTABLE' ? 'warn' : 'bad'
  return <span className={`pill ${kind}`}>{value.toLowerCase()}</span>
}
