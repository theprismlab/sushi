import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import { useNavigate, useSearchParams } from 'react-router-dom'
import { api } from './api.js'
import { useCatalog, useIdentity } from './App.jsx'

const MODULE_GROUPS = new Set(['modules', 'qc_modules', 'analytics', 'portal', 'metadata'])

export default function Launch() {
  const catalog = useCatalog()
  const navigate = useNavigate()
  const [searchParams] = useSearchParams()
  const [who] = useIdentity()

  const fromRun = searchParams.get('from')
  const [preset, setPreset] = useState(searchParams.get('preset') || catalog.presets[0].id)
  // Params of a previous run being re-run. Applied once, on top of the preset.
  const [pending, setPending] = useState(null)
  const [values, setValues] = useState(null)
  const [fromPreset, setFromPreset] = useState([])
  const [buildDirs, setBuildDirs] = useState(null)
  const [touched, setTouched] = useState(() => new Set())
  const [showAdvanced, setShowAdvanced] = useState(false)
  const [check, setCheck] = useState(null)
  const [checking, setChecking] = useState(false)
  const [archive, setArchive] = useState(false)
  const [launching, setLaunching] = useState(false)
  const [launchError, setLaunchError] = useState(null)

  const stock = useMemo(
    () => Object.fromEntries(catalog.params.map((p) => [p.name, p.default])),
    [catalog],
  )
  const byGroup = useMemo(() => {
    const map = new Map(catalog.groups.map((g) => [g.id, { ...g, params: [] }]))
    catalog.params.forEach((p) => map.get(p.group).params.push(p))
    return [...map.values()]
  }, [catalog])

  useEffect(() => {
    if (!fromRun) return
    api.run(fromRun).then((r) => { setPreset(r.preset); setPending(r.params) })
  }, [fromRun])

  // Selecting a screen type is the one action that repopulates the whole form.
  useEffect(() => {
    if (fromRun && !pending) return // wait for the run we are copying
    let cancelled = false
    setValues(null)
    setTouched(new Set())
    setArchive(false)
    Promise.all([api.presetDefaults(preset, ''), api.buildDirs(preset)]).then(([d, dirs]) => {
      if (cancelled) return
      setValues(pending ? { ...d.values, ...pending } : d.values)
      setFromPreset(d.from_preset)
      setBuildDirs(dirs)
    })
    return () => { cancelled = true }
  }, [preset, pending, fromRun])

  // Clicking a screen type is also how you discard a copied run's overrides.
  function choosePreset(id) {
    setPending(null)
    setPreset(id)
  }

  const setValue = useCallback((name, value) => {
    setTouched((t) => new Set(t).add(name))
    setValues((v) => {
      const next = { ...v, [name]: value }
      // Build name drives screen and build dir until the user overrides them.
      if (name === 'BUILD_NAME') {
        const dir = catalog.presets.find((p) => p.id === preset).dir
        if (!touched.has('SCREEN')) next.SCREEN = value
        if (!touched.has('BUILD_DIR')) next.BUILD_DIR = `${catalog.prismseq_root}/${dir}/${value}`
      }
      return next
    })
  }, [catalog, preset, touched])

  // Preflight is cheap (stat calls and a json read) so run it as they type.
  const request = useMemo(
    () => (values ? { preset, launched_by: who || 'unknown', values, archive_existing_config: archive } : null),
    [preset, who, values, archive],
  )
  const latest = useRef(0)
  useEffect(() => {
    if (!request) return
    const token = ++latest.current
    setChecking(true)
    const timer = setTimeout(() => {
      api.preflight(request)
        .then((r) => { if (token === latest.current) setCheck(r) })
        .catch(() => { if (token === latest.current) setCheck(null) })
        .finally(() => { if (token === latest.current) setChecking(false) })
    }, 400)
    return () => clearTimeout(timer)
  }, [request])

  async function launch() {
    setLaunching(true)
    setLaunchError(null)
    try {
      const { id } = await api.launch(request)
      navigate(`/runs/${id}`)
    } catch (e) {
      setLaunchError(e.detail?.problems?.join(' ') || e.message)
      setLaunching(false)
    }
  }

  if (!values) return <div className="muted">Loading defaults…</div>

  const selected = catalog.presets.find((p) => p.id === preset)
  const blocked = !who.trim() || (check?.problems?.length ?? 0) > 0
  const dirEntry = buildDirs?.builds?.find((b) => b.name === values.BUILD_NAME)

  return (
    <div className="launch">
      <section className="card">
        <h2>1. Screen type</h2>
        <p className="muted">
          Sets every parameter below. Defaults come from the most recent build of each type.
          {pending && ' Currently showing a copy of an earlier run — click a screen type to discard those overrides.'}
        </p>
        <div className="preset-grid">
          {catalog.presets.map((p) => (
            <button
              key={p.id}
              className={`preset ${p.id === preset ? 'on' : ''}`}
              onClick={() => choosePreset(p.id)}
            >
              <strong>{p.id}</strong>
              <span className="muted">{p.label.replace(/^[A-Z]+ — /, '')}</span>
              {p.reference && <span className="tag">from {p.reference}</span>}
            </button>
          ))}
        </div>
      </section>

      <section className="card">
        <h2>2. Build</h2>
        <div className="field">
          <label htmlFor="build-name">
            Build name <span className="req">required</span>
          </label>
          <input
            id="build-name"
            list="existing-builds"
            value={values.BUILD_NAME}
            onChange={(e) => setValue('BUILD_NAME', e.target.value)}
            placeholder={selected.reference ? `e.g. ${selected.reference}` : 'e.g. MTS033'}
            autoComplete="off"
          />
          <datalist id="existing-builds">
            {buildDirs?.builds?.map((b) => <option key={b.name} value={b.name} />)}
          </datalist>
          <p className="help">
            {buildDirs && !buildDirs.exists && <>Directory {buildDirs.root} does not exist. </>}
            {dirEntry ? (
              <>
                Existing build, last modified {dirEntry.modified.slice(0, 10)}.{' '}
                {dirEntry.has_raw_counts ? 'Raw counts present.' : 'No raw counts file.'}
                {dirEntry.has_config && ' Has a config.json.'}
              </>
            ) : (
              <>{buildDirs?.builds?.length ?? 0} existing builds under {buildDirs?.root}</>
            )}
          </p>
        </div>

        {byGroup
          .filter((g) => g.id === 'build')
          .map((g) => (
            <div className="grid" key={g.id}>
              {g.params
                .filter((p) => p.name !== 'BUILD_NAME')
                .map((p) => (
                  <Field key={p.name} spec={p} value={values[p.name]} stock={stock[p.name]}
                         fromPreset={fromPreset.includes(p.name)} onChange={setValue} />
                ))}
            </div>
          ))}
      </section>

      <section className="card">
        <h2>3. Modules</h2>
        <p className="muted">
          Which scripts run. Jenkins does not record these in <code>config.json</code>, so this UI is
          the only place the choice is kept.
        </p>
        {byGroup
          .filter((g) => MODULE_GROUPS.has(g.id))
          .map((g) => (
            <div className="module-block" key={g.id}>
              <h3>{g.label}</h3>
              {g.help && <p className="help">{g.help}</p>}
              <div className="module-grid">
                {g.params.map((p) => (
                  <Toggle key={p.name} spec={p} value={values[p.name]} stock={stock[p.name]}
                          fromPreset={fromPreset.includes(p.name)} onChange={setValue} />
                ))}
              </div>
            </div>
          ))}
      </section>

      <section className="card">
        <h2>4. Controls & signatures</h2>
        {byGroup
          .filter((g) => g.id === 'controls')
          .map((g) => (
            <div className="grid" key={g.id}>
              {g.params.map((p) => (
                <Field key={p.name} spec={p} value={values[p.name]} stock={stock[p.name]}
                       fromPreset={fromPreset.includes(p.name)} onChange={setValue} />
              ))}
            </div>
          ))}
      </section>

      <section className="card">
        <button className="disclosure" onClick={() => setShowAdvanced((s) => !s)}>
          {showAdvanced ? '▾' : '▸'} Advanced parameters
          <span className="muted">
            {' '}
            {byGroup.filter((g) => g.tier === 'advanced').reduce((n, g) => n + g.params.length, 0)}{' '}
            settings — filenames, column names, QC thresholds, pipeline version
          </span>
        </button>
        {showAdvanced &&
          byGroup
            .filter((g) => g.tier === 'advanced')
            .map((g) => (
              <div className="module-block" key={g.id}>
                <h3>{g.label}</h3>
                {g.help && <p className="help">{g.help}</p>}
                <div className="grid">
                  {g.params.map((p) => (
                    <Field key={p.name} spec={p} value={values[p.name]} stock={stock[p.name]}
                           fromPreset={fromPreset.includes(p.name)} onChange={setValue} />
                  ))}
                </div>
              </div>
            ))}
      </section>

      <section className="card sticky">
        <h2>Launch</h2>
        {checking && <p className="muted">Checking…</p>}

        {check?.problems?.length > 0 && (
          <div className="banner error">
            <strong>Fix before launching</strong>
            <ul>{check.problems.map((p) => <li key={p}>{p}</li>)}</ul>
          </div>
        )}

        {check?.stale_config && (
          <div className="banner warn">
            <strong>An existing config.json will override the form</strong>
            <p>{check.stale_config.note}</p>
            <table className="mini">
              <thead><tr><th>File</th><th>Parameter</th><th>On disk (wins)</th><th>Your value</th></tr></thead>
              <tbody>
                {check.stale_config.conflicts.map((c) => (
                  <tr key={`${c.file}:${c.name}`}>
                    <td><code>{c.file}</code></td>
                    <td><code>{c.name}</code></td>
                    <td>{c.existing}</td>
                    <td>{c.submitted}</td>
                  </tr>
                ))}
              </tbody>
            </table>
            <label className="checkbox">
              <input type="checkbox" checked={archive} onChange={(e) => setArchive(e.target.checked)} />
              Archive the existing config.json and qc_params.json (copied to <code>.bak</code>) so
              these values take effect
            </label>
          </div>
        )}

        {!who.trim() && <p className="banner warn">Enter your name in the header first — it is recorded against the run.</p>}
        {launchError && <div className="banner error">{launchError}</div>}

        <div className="launch-row">
          <button className="primary" disabled={blocked || launching} onClick={launch}>
            {launching ? 'Triggering…' : `Launch ${values.BUILD_NAME || 'build'}`}
          </button>
          <span className="muted">
            {values.BUILD_DIR} · branch {values.GIT_BRANCH}
            {!values.TRIGGER_BUILD && ' · config only, no modules will run'}
          </span>
        </div>
      </section>
    </div>
  )
}

function changedClass(value, stock, fromPreset) {
  if (String(value) !== String(stock)) return fromPreset ? 'preset-set' : 'user-set'
  return ''
}

function Field({ spec, value, stock, fromPreset, onChange }) {
  const id = `f-${spec.name}`
  const mark = changedClass(value, stock, fromPreset)
  return (
    <div className={`field ${mark}`}>
      <label htmlFor={id}>
        {spec.label || spec.name}
        {spec.required && <span className="req">required</span>}
        {mark === 'preset-set' && <span className="tag">screen default</span>}
        {mark === 'user-set' && <span className="tag alt">edited</span>}
      </label>
      {spec.type === 'choice' ? (
        <select id={id} value={value} onChange={(e) => onChange(spec.name, e.target.value)}>
          {spec.choices.map((c) => <option key={c} value={c}>{c}</option>)}
        </select>
      ) : spec.type === 'bool' ? (
        <input id={id} type="checkbox" checked={!!value}
               onChange={(e) => onChange(spec.name, e.target.checked)} />
      ) : (
        <input id={id} value={value ?? ''} onChange={(e) => onChange(spec.name, e.target.value)} />
      )}
      <code className="pname">{spec.name}</code>
      {spec.help && <p className="help">{spec.help}</p>}
    </div>
  )
}

function Toggle({ spec, value, stock, fromPreset, onChange }) {
  const mark = changedClass(value, stock, fromPreset)
  return (
    <label className={`toggle ${value ? 'on' : ''} ${mark}`} title={spec.help || spec.name}>
      <input type="checkbox" checked={!!value} onChange={(e) => onChange(spec.name, e.target.checked)} />
      <span>
        {spec.label}
        {mark === 'preset-set' && <span className="tag">screen default</span>}
      </span>
      <code className="pname">{spec.name}</code>
    </label>
  )
}
