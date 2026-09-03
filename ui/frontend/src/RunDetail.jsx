import { useEffect, useMemo, useRef, useState } from 'react'
import { Link, useParams } from 'react-router-dom'
import { api, formatBytes, formatDuration, formatTime, RUNNING_STATES } from './api.js'
import { useCatalog } from './App.jsx'
import { Status } from './History.jsx'

export default function RunDetail() {
  const { id } = useParams()
  const catalog = useCatalog()
  const [run, setRun] = useState(null)
  const [outputs, setOutputs] = useState(null)
  const [error, setError] = useState(null)

  const live = run && RUNNING_STATES.includes(run.status)

  useEffect(() => {
    let cancelled = false
    const load = () => api.run(id)
      .then((r) => { if (!cancelled) { setRun(r); setError(null) } })
      .catch((e) => { if (!cancelled) setError(e.message) })
    load()
    const timer = setInterval(load, 5000)
    return () => { cancelled = true; clearInterval(timer) }
  }, [id])

  useEffect(() => {
    if (!run) return
    api.outputs(id).then(setOutputs).catch(() => setOutputs(null))
  }, [id, run?.status])

  if (error) return <div className="banner error">{error}</div>
  if (!run) return <div className="muted">Loading run…</div>

  return (
    <div className="detail">
      <section className="card head">
        <div>
          <h2>
            {run.build_name} <Status value={run.status} />
          </h2>
          <p className="muted">
            {run.preset}
            {run.screen_type !== run.preset && ` (sent as SCREEN_TYPE=${run.screen_type})`}
            {' · '}launched by {run.launched_by} · {formatTime(run.created_at)}
          </p>
          <p className="muted small"><code>{run.build_dir}</code></p>
        </div>
        <div className="head-actions">
          {run.jenkins_url && (
            <a className="button" href={run.jenkins_url} target="_blank" rel="noreferrer">Open in Jenkins</a>
          )}
          <Link className="button" to={`/launch?from=${run.id}`}>Re-run with these settings</Link>
          {live && run.build_number && (
            <button onClick={() => api.stop(run.id).then(setRun).catch((e) => setError(e.message))}>
              Stop build
            </button>
          )}
        </div>
      </section>

      {run.error && <div className="banner error">{run.error}</div>}

      <div className="stat-row card">
        <Stat label="Jenkins build" value={run.build_number ? `#${run.build_number}` : 'queued'} />
        <Stat label="Started" value={formatTime(run.started_at)} />
        <Stat label="Finished" value={formatTime(run.finished_at)} />
        <Stat label="Duration" value={formatDuration(run.duration_ms)} />
        <Stat label="Pipeline branch" value={run.git_branch || '—'} />
        <Stat label="Modules run" value={run.enabled_modules.length} />
      </div>

      <Notes run={run} onSaved={setRun} />

      <section className="card">
        <h3>Modules</h3>
        <div className="module-grid readonly">
          {Object.entries(run.modules).sort().map(([name, on]) => (
            <span key={name} className={`chip ${on ? 'on' : 'off'}`}>
              {on ? '✓' : '·'} {name}
            </span>
          ))}
        </div>
      </section>

      <Outputs outputs={outputs} />
      <Params run={run} catalog={catalog} />
      <Log runId={id} live={live} status={run.status} />
    </div>
  )
}

function Stat({ label, value }) {
  return <div className="stat"><span className="muted">{label}</span><strong>{value}</strong></div>
}

function Notes({ run, onSaved }) {
  const [text, setText] = useState(run.notes || '')
  const [saved, setSaved] = useState(false)
  useEffect(() => { setText(run.notes || '') }, [run.id])
  return (
    <section className="card">
      <h3>Notes</h3>
      <textarea
        rows={3}
        value={text}
        placeholder="What was this run for, what came out of it, anything odd."
        onChange={(e) => { setText(e.target.value); setSaved(false) }}
      />
      <div className="launch-row">
        <button onClick={() => api.setNotes(run.id, text).then((r) => { onSaved(r); setSaved(true) })}>
          Save notes
        </button>
        {saved && <span className="muted">Saved.</span>}
      </div>
    </section>
  )
}

function Outputs({ outputs }) {
  if (!outputs) return null
  return (
    <section className="card">
      <h3>Outputs</h3>
      {!outputs.exists ? (
        <p className="muted">Build directory {outputs.build_dir} is not readable from this host.</p>
      ) : (
        <>
          <table className="mini full">
            <thead><tr><th>File</th><th>Size</th><th>Modified</th></tr></thead>
            <tbody>
              {outputs.files.map((f) => (
                <tr key={f.param} className={f.exists ? '' : 'absent'}>
                  <td><code>{f.name}</code></td>
                  <td>{f.exists ? formatBytes(f.size) : <span className="muted">not written</span>}</td>
                  <td className="small">
                    {f.modified ? f.modified.replace('T', ' ').slice(0, 16) : '—'}
                    {f.stale && <span className="tag alt" title="Older than this run — left over from a previous build">stale</span>}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
          {outputs.qc_tables.length > 0 && (
            <p className="help">qc_tables/: {outputs.qc_tables.join(', ')}</p>
          )}
        </>
      )}
    </section>
  )
}

function Params({ run, catalog }) {
  const [showAll, setShowAll] = useState(false)
  const index = useMemo(() => Object.fromEntries(catalog.params.map((p) => [p.name, p])), [catalog])
  const rows = Object.entries(run.params)
    .filter(([name]) => index[name] && index[name].type !== 'bool')
    .map(([name, value]) => ({ name, value, spec: index[name], changed: String(value) !== String(index[name].default) }))
    .sort((a, b) => Number(b.changed) - Number(a.changed) || a.name.localeCompare(b.name))
  const visible = showAll ? rows : rows.filter((r) => r.changed)

  return (
    <section className="card">
      <h3>Parameters</h3>
      <p className="help">
        Non-default values first. Booleans are shown under Modules above.
      </p>
      <table className="mini full">
        <thead><tr><th>Parameter</th><th>Value</th><th>Default</th></tr></thead>
        <tbody>
          {visible.map((r) => (
            <tr key={r.name} className={r.changed ? 'changed' : ''}>
              <td><code>{r.name}</code></td>
              <td className="wrap">{String(r.value) || <span className="muted">empty</span>}</td>
              <td className="muted wrap small">{r.changed ? String(r.spec.default) || 'empty' : ''}</td>
            </tr>
          ))}
        </tbody>
      </table>
      <button className="disclosure" onClick={() => setShowAll((s) => !s)}>
        {showAll ? 'Show only non-default' : `Show all ${rows.length} parameters`}
      </button>
    </section>
  )
}

function Log({ runId, live, status }) {
  const [text, setText] = useState('')
  const [next, setNext] = useState(0)
  const [message, setMessage] = useState(null)
  const [follow, setFollow] = useState(true)
  const box = useRef(null)
  const offset = useRef(0)

  useEffect(() => { setText(''); setNext(0); offset.current = 0 }, [runId])

  useEffect(() => {
    let cancelled = false
    const poll = () => api.log(runId, offset.current)
      .then((r) => {
        if (cancelled) return
        if (r.text) {
          setText((t) => t + r.text)
          offset.current = r.next
          setNext(r.next)
        }
        setMessage(r.message || null)
      })
      .catch(() => {})
    poll()
    if (!live) return () => { cancelled = true }
    const timer = setInterval(poll, 3000)
    return () => { cancelled = true; clearInterval(timer) }
  }, [runId, live, status])

  useEffect(() => {
    if (follow && box.current) box.current.scrollTop = box.current.scrollHeight
  }, [text, follow])

  return (
    <section className="card">
      <h3>
        Console log
        {live && <span className="pill busy">streaming</span>}
      </h3>
      {message && <p className="help">{message}</p>}
      <label className="checkbox">
        <input type="checkbox" checked={follow} onChange={(e) => setFollow(e.target.checked)} />
        Follow output
      </label>
      <pre className="log" ref={box}>{text || (live ? 'Waiting for output…' : 'No log available.')}</pre>
      <p className="muted small">{next.toLocaleString()} bytes</p>
    </section>
  )
}
