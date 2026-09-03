import { createContext, useCallback, useContext, useEffect, useMemo, useRef, useState } from 'react'
import { useNavigate, useSearchParams } from 'react-router-dom'
import { api } from './api.js'
import { useCatalog, useIdentity } from './App.jsx'

const MODULE_GROUPS = new Set(['modules', 'qc_modules', 'analytics', 'portal', 'metadata'])

// cellDB-backed value lists, keyed by param name. Via context so every Field
// picks them up without threading a prop through four call sites.
const SuggestContext = createContext(null)

// BUILD_DIR is owned by the path picker in step 1, so it is not rendered as a
// free-text field; BUILD_NAME has its own input above the rest of the group.
const BUILD_FIELDS_HANDLED_ELSEWHERE = new Set(['BUILD_NAME', 'BUILD_DIR'])

/** Strip PRISMSEQ_ROOT off an absolute build dir to get the relative path. */
function relativeToRoot(dir, root) {
  if (!dir) return ''
  const prefix = `${root.replace(/\/+$/, '')}/`
  return dir.startsWith(prefix) ? dir.slice(prefix.length).replace(/\/+$/, '') : ''
}

export default function Launch() {
  const catalog = useCatalog()
  const navigate = useNavigate()
  const [searchParams] = useSearchParams()
  const [identity] = useIdentity()

  const root = catalog.prismseq_root
  const fromRun = searchParams.get('from')

  // The build directory is now the first thing entered: the screen type and
  // the name/screen defaults are all derived from it.
  const [buildPath, setBuildPath] = useState('')
  // Set only when the user overrides the screen type inferred from the path.
  const [presetOverride, setPresetOverride] = useState(searchParams.get('preset') || '')
  const [resolved, setResolved] = useState(null)
  const [paths, setPaths] = useState(null)
  // Real values from cellDB for the few params that have them.
  const [suggestions, setSuggestions] = useState(null)

  // Params of a previous run being re-run. Applied once, on top of the preset.
  const [pending, setPending] = useState(null)
  const [values, setValues] = useState(null)
  const [fromPreset, setFromPreset] = useState([])
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

  useEffect(() => { api.buildPaths().then(setPaths).catch(() => setPaths(null)) }, [])
  useEffect(() => { api.suggestions().then(setSuggestions).catch(() => setSuggestions({})) }, [])

  useEffect(() => {
    if (!fromRun) return
    api.run(fromRun).then((r) => {
      const { BUILD_DIR, ...rest } = r.params
      // The path picker owns BUILD_DIR; keep one source of truth for it.
      setBuildPath(relativeToRoot(BUILD_DIR, root))
      setPresetOverride(r.preset)
      setPending(rest)
    })
  }, [fromRun, root])

  // One resolve for both inputs: the path infers a screen type, an explicit
  // screen type overrides that, and either repopulates the whole form.
  useEffect(() => {
    if (fromRun && !pending) return // wait for the run we are copying
    let cancelled = false
    const timer = setTimeout(() => {
      api.pathDefaults(buildPath, presetOverride).then((d) => {
        if (cancelled) return
        setResolved(d)
        setFromPreset(d.from_preset)
        setValues(pending ? { ...d.values, ...pending } : d.values)
        setTouched(new Set())
        setArchive(false)
      }).catch(() => { if (!cancelled) setResolved(null) })
    }, 250)
    return () => { cancelled = true; clearTimeout(timer) }
  }, [buildPath, presetOverride, pending, fromRun])

  // Typing a new path re-infers the screen type; an earlier manual override
  // would otherwise stick to a directory it no longer matches.
  function choosePath(next) {
    setBuildPath(next)
    setPresetOverride('')
    setPending(null)
  }

  // Clicking a screen type is also how you discard a copied run's overrides.
  function choosePreset(id) {
    setPending(null)
    setPresetOverride(id)
  }

  const setValue = useCallback((name, value) => {
    setTouched((t) => new Set(t).add(name))
    setValues((v) => {
      const next = { ...v, [name]: value }
      // Build name still seeds SCREEN, but no longer rewrites BUILD_DIR --
      // the path picker owns that now.
      if (name === 'BUILD_NAME' && !touched.has('SCREEN')) next.SCREEN = value
      return next
    })
  }, [touched])

  // Preflight is cheap (stat calls and a json read) so run it as they type --
  // but not before a build directory is picked. Until then BUILD_DIR is still
  // PRISMSEQ_ROOT, so every check is about the wrong directory: the required
  // fields are empty because nothing is selected, and any stray config.json at
  // the root reports a wall of overrides that has nothing to do with the build.
  const selected = Boolean(buildPath.trim())
  const preset = resolved?.preset
  const request = useMemo(
    () => (values && preset && selected
      ? { preset, values, archive_existing_config: archive }
      : null),
    [preset, values, archive, selected],
  )
  const latest = useRef(0)
  useEffect(() => {
    if (!request) {
      // Clearing the path must also clear the last result, or the warnings
      // from the previous selection linger over an empty form.
      latest.current += 1
      setCheck(null)
      setChecking(false)
      return
    }
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

  if (!values || !resolved) return <div className="muted">Loading defaults…</div>

  const blocked = !identity || !buildPath.trim() || (check?.problems?.length ?? 0) > 0
  const entry = paths?.paths?.find((p) => p.path === resolved.path)

  return (
    <SuggestContext.Provider value={suggestions}>
    <div className="launch">
      <section className="card">
        <h2>1. Build directory</h2>
        <p className="muted">
          Browse to the build directory. Everything else — screen type, build name, screen — is
          derived from it.
        </p>
        <BuildPicker
          root={root}
          value={resolved.path}
          allPaths={paths?.paths}
          onSelect={choosePath}
        />
        <p className="help">
          {!buildPath.trim() ? (
            <>Nothing selected yet. {paths?.paths?.length ?? 0} build directories exist.</>
          ) : resolved.exists ? (
            <>
              <code>{resolved.build_dir}</code>
              {entry && <> · last modified {entry.modified.slice(0, 10)}</>}
              {' · '}
              {resolved.has_raw_counts ? 'raw counts present' : 'no raw counts file'}
              {resolved.has_config && ' · has a config.json'}
            </>
          ) : (
            <>
              <span className="tag alt">missing</span> <code>{resolved.build_dir}</code> does not
              exist. The pipeline needs the directory and its input files to be there already.
            </>
          )}
        </p>
      </section>

      <section className="card">
        <h2>
          2. Screen type
          {resolved.inferred && <span className="tag">from the path</span>}
        </h2>
        <p className="muted">
          {resolved.inferred
            ? `Inferred from ${resolved.preset_reason}. Click another to override.`
            : `Set manually. Clear the path or retype it to go back to inferring.`}
          {pending && ' Currently showing a copy of an earlier run — click a screen type to discard those overrides.'}
        </p>
        <div className="preset-grid">
          {catalog.presets.map((p) => (
            <button
              key={p.id}
              className={`preset ${p.id === resolved.preset ? 'on' : ''}`}
              onClick={() => choosePreset(p.id)}
            >
              <strong>{p.id}</strong>
              <span className="muted">{p.label.replace(/^[A-Z]+ — /, '')}</span>
              {/* Deliberately not p.reference: defaults come from the screen
                  types page, not from that build. */}
              <span className="muted">
                {p.override_count} default{p.override_count === 1 ? '' : 's'}
              </span>
            </button>
          ))}
        </div>
      </section>

      <section className="card">
        <h2>3. Build</h2>
        <div className="field">
          <label htmlFor="build-name">
            Build name <span className="req">required</span>
          </label>
          <input
            id="build-name"
            value={values.BUILD_NAME}
            onChange={(e) => setValue('BUILD_NAME', e.target.value)}
            placeholder="e.g. MTS033"
            autoComplete="off"
          />
          <p className="help">
            Defaults to the last segment of the build directory. Passed to Jenkins as{' '}
            <code>BUILD_NAME</code>.
          </p>
        </div>

        {byGroup
          .filter((g) => g.id === 'build')
          .map((g) => (
            <div className="grid" key={g.id}>
              {g.params
                .filter((p) => !BUILD_FIELDS_HANDLED_ELSEWHERE.has(p.name))
                .map((p) => (
                  <Field key={p.name} spec={p} value={values[p.name]} stock={stock[p.name]}
                         fromPreset={fromPreset.includes(p.name)} onChange={setValue} />
                ))}
            </div>
          ))}
      </section>

      <section className="card">
        <h2>4. Modules</h2>
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
        <h2>5. Controls & signatures</h2>
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
        <h2>6. Pipeline version</h2>
        <PipelineVersion values={values} stock={stock} fromPreset={fromPreset} onChange={setValue} />
      </section>

      <section className="card">
        <button className="disclosure" onClick={() => setShowAdvanced((s) => !s)}>
          {showAdvanced ? '▾' : '▸'} Advanced parameters
          <span className="muted">
            {' '}
            {byGroup.filter((g) => g.tier === 'advanced').reduce((n, g) => n + g.params.length, 0)}{' '}
            settings — filenames, column names, QC thresholds
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
        {!selected && (
          <p className="muted">
            Pick a build directory in step 1 first. Nothing is checked until then.
          </p>
        )}
        {selected && checking && <p className="muted">Checking…</p>}

        {check?.problems?.length > 0 && (
          <div className="banner error">
            <strong>Fix before launching</strong>
            <ul>{check.problems.map((p) => <li key={p}>{p}</li>)}</ul>
          </div>
        )}

        {check?.warnings?.length > 0 && (
          <div className="banner warn">
            <strong>Worth knowing</strong>
            <ul>{check.warnings.map((w) => <li key={w}>{w}</li>)}</ul>
          </div>
        )}

        {check?.stale_config && (
          <div className="banner warn">
            <strong>An existing config.json will override the form</strong>
            <p>{check.stale_config.note}</p>
            <div className="table-scroll">
              <table className="mini conflicts">
                <thead><tr><th>File</th><th>Parameter</th><th>On disk (wins)</th><th>Your value</th></tr></thead>
                <tbody>
                  {check.stale_config.conflicts.map((c) => (
                    <tr key={`${c.file}:${c.name}`}>
                      <td><code className="wrap">{c.file}</code></td>
                      <td><code className="wrap">{c.name}</code></td>
                      <td className="wrap">{c.existing}</td>
                      <td className="wrap">{c.submitted}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            <label className="checkbox">
              <input type="checkbox" checked={archive} onChange={(e) => setArchive(e.target.checked)} />
              Archive the existing config.json and qc_params.json (copied to <code>.bak</code>) so
              these values take effect
            </label>
          </div>
        )}

        {!identity && (
          <p className="banner warn">
            Sign in to Jenkins from the header first. Jenkins requires it to queue a build, and it
            is what the run log records as the launcher.
          </p>
        )}
        {launchError && <div className="banner error">{launchError}</div>}

        <div className="launch-row">
          <button className="primary" disabled={blocked || launching} onClick={launch}>
            {launching ? 'Triggering…' : `Launch ${values.BUILD_NAME || 'build'}`}
          </button>
          <span className="muted wrap">
            {selected ? (
              <>
                {values.BUILD_DIR} · branch {values.GIT_BRANCH}
                {!values.TRIGGER_BUILD && ' · config only, no modules will run'}
              </>
            ) : (
              // Showing PRISMSEQ_ROOT here reads as if the root were the target.
              'No build directory selected'
            )}
          </span>
        </div>
      </section>
    </div>
    </SuggestContext.Provider>
  )
}

/** Directory browser over the real prismSeq listings.
 *
 *  Keyboard first: the filter narrows the current level as you type, Enter
 *  opens a folder (or selects a build), and once you are inside a folder the
 *  filter clears so you see everything in it. Typing two or more characters
 *  also searches every build path, so you can jump straight to a build name
 *  without knowing which screen-type folder it lives under.
 *
 *  A build directory is a leaf you select; anything else is a folder you open.
 *  Some build directories also contain subdirectories, so those offer both. */
function BuildPicker({ root, value, allPaths, onSelect }) {
  // Open at the folder holding the current selection, not at the root.
  const [cwd, setCwd] = useState(() => (value ? parentPath(value) : ''))
  const [listing, setListing] = useState(null)
  const [loading, setLoading] = useState(true)
  const [query, setQuery] = useState('')
  const [cursor, setCursor] = useState(0)

  useEffect(() => {
    let cancelled = false
    setLoading(true)
    api.ls(cwd)
      .then((d) => { if (!cancelled) setListing(d) })
      .catch(() => { if (!cancelled) setListing(null) })
      .finally(() => { if (!cancelled) setLoading(false) })
    return () => { cancelled = true }
  }, [cwd])

  const q = query.trim().toLowerCase()

  // Entries at this level whose name matches, plus — for a query long enough
  // to be meaningful — builds anywhere else in the tree.
  const here = useMemo(
    () => (listing?.entries || []).filter((e) => !q || e.name.toLowerCase().includes(q)),
    [listing, q],
  )
  const elsewhere = useMemo(() => {
    if (q.length < 2 || !allPaths) return []
    const shown = new Set(here.map((e) => e.path))
    return allPaths
      .filter((p) => !shown.has(p.path) && p.path.toLowerCase().includes(q))
      .slice(0, 40)
  }, [allPaths, q, here])

  const items = useMemo(
    () => [
      ...here.map((e) => ({ kind: e.is_build ? 'build' : 'dir', path: e.path, label: e.name, entry: e })),
      ...elsewhere.map((p) => ({ kind: 'build', path: p.path, label: p.path, entry: p, remote: true })),
    ],
    [here, elsewhere],
  )

  // Any change to what is on screen invalidates the highlight position.
  useEffect(() => { setCursor(0) }, [cwd, q])

  function activate(item) {
    if (!item) return
    if (item.kind === 'dir') {
      setCwd(item.path)
      setQuery('') // now show everything inside, which is the point of opening it
    } else {
      onSelect(item.path)
    }
  }

  function onKeyDown(e) {
    if (e.key === 'ArrowDown' || (e.key === 'Tab' && !e.shiftKey && items.length)) {
      e.preventDefault()
      setCursor((c) => Math.min(c + 1, items.length - 1))
    } else if (e.key === 'ArrowUp' || (e.key === 'Tab' && e.shiftKey && items.length)) {
      e.preventDefault()
      setCursor((c) => Math.max(c - 1, 0))
    } else if (e.key === 'Enter') {
      e.preventDefault()
      activate(items[cursor])
    } else if (e.key === 'Escape') {
      e.preventDefault()
      setQuery('')
    } else if (e.key === 'Backspace' && !query && listing?.parent !== null) {
      e.preventDefault()
      setCwd(listing.parent) // empty filter + backspace reads as "go up"
    }
  }

  const segments = cwd ? cwd.split('/') : []

  return (
    <div className="browser">
      <div className="crumbs">
        <button className="crumb" onClick={() => { setCwd(''); setQuery('') }} title={root}>
          prismSeq
        </button>
        {segments.map((seg, i) => (
          <span key={i}>
            <span className="crumb-sep">/</span>
            <button className="crumb"
                    onClick={() => { setCwd(segments.slice(0, i + 1).join('/')); setQuery('') }}>
              {seg}
            </button>
          </span>
        ))}
      </div>

      <div className="picker-filter">
        <input
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          onKeyDown={onKeyDown}
          placeholder="Type to filter, ↑↓ to move, Enter to open or select"
          aria-label="Filter build directories"
          autoComplete="off"
          spellCheck="false"
        />
        <span className="muted small">
          {loading ? 'reading…' : `${here.length}${q ? ` of ${listing?.entries?.length ?? 0}` : ''} here`}
          {elsewhere.length > 0 && ` · ${elsewhere.length} elsewhere`}
        </span>
      </div>

      <div className="entries">
        {loading && <p className="muted entry-msg">Reading {cwd || 'prismSeq'}…</p>}

        {!loading && listing?.parent !== null && !q && (
          <div className="entry-row">
            <button className="entry" onClick={() => setCwd(listing.parent)}>
              <span className="entry-icon">↑</span>
              <span className="entry-name">..</span>
            </button>
          </div>
        )}

        {!loading && items.length === 0 && (
          <p className="muted entry-msg">
            {q ? `Nothing matching “${query}” here or anywhere else.` : 'No subdirectories here.'}
          </p>
        )}

        {!loading && items.map((item, i) => (
          <div
            key={`${item.remote ? 'r' : 'l'}:${item.path}`}
            className={`entry-row ${item.path === value ? 'selected' : ''} ${i === cursor ? 'active' : ''}`}
            ref={i === cursor ? scrollIntoView : null}
          >
            <button
              className="entry"
              onMouseEnter={() => setCursor(i)}
              onClick={() => activate(item)}
              title={item.kind === 'build' ? `Select ${item.path}` : `Open ${item.path}`}
            >
              <span className="entry-icon">{item.kind === 'build' ? '▣' : '▸'}</span>
              <span className="entry-name">{item.label}</span>
              {item.kind === 'build' && <span className="tag">build</span>}
              {item.kind === 'build' && !item.entry.has_raw_counts && (
                <span className="tag alt">no counts</span>
              )}
              <span className="entry-date muted">{item.entry.modified.slice(0, 10)}</span>
            </button>
            {/* A build directory can still hold subdirectories worth opening. */}
            {item.kind === 'build' && item.entry.child_dirs > 0 && (
              <button className="entry-open" onClick={() => { setCwd(item.path); setQuery('') }}
                      title={`Open ${item.label}`}>
                open
              </button>
            )}
          </div>
        ))}
      </div>
    </div>
  )
}

function scrollIntoView(el) {
  el?.scrollIntoView({ block: 'nearest' })
}

function parentPath(rel) {
  const parts = (rel || '').split('/')
  return parts.slice(0, -1).join('/')
}

function changedClass(value, stock, fromPreset) {
  if (String(value) !== String(stock)) return fromPreset ? 'preset-set' : 'user-set'
  return ''
}

/** GIT_BRANCH / USE_LATEST / COMMIT_ID, populated from GitHub.
 *
 *  These three decide which code the build runs, which is why they get their
 *  own section instead of sitting behind the advanced disclosure. Typed
 *  freehand, a typo is only found when Jenkins fails to check out -- or, worse,
 *  when it silently runs a revision nobody meant.
 *
 *  Both fields stay editable text if GitHub cannot be reached, so an outage
 *  degrades the section rather than blocking a launch.
 */
function PipelineVersion({ values, stock, fromPreset, onChange }) {
  const branch = values.GIT_BRANCH || ''
  const useLatest = !!values.USE_LATEST
  const [refs, setRefs] = useState(null)

  useEffect(() => {
    let cancelled = false
    // Commits are fetched for the selected branch, so this refetches on change.
    // The backend caches, so flipping between branches is not 20 API calls.
    api.gitRefs(branch)
      .then((r) => { if (!cancelled) setRefs(r) })
      .catch(() => { if (!cancelled) setRefs({ branches: [], commits: [], error: 'lookup failed' }) })
    return () => { cancelled = true }
  }, [branch])

  const branches = refs?.branches || []
  const commits = refs?.commits || []
  const head = commits[0]
  // A branch set from a config.json or an old run may not exist any more.
  // Keep it selectable and say so rather than silently snapping to another.
  const unknownBranch = Boolean(branch && branches.length && !branches.includes(branch))

  return (
    <>
      {refs?.error && (
        <p className="help warn-text">
          Could not read branches from GitHub ({refs.error}). Both fields are still editable.
        </p>
      )}

      <div className="grid">
        <div className={`field ${changedClass(branch, stock.GIT_BRANCH, fromPreset.includes('GIT_BRANCH'))}`}>
          <label htmlFor="f-GIT_BRANCH">
            Pipeline branch
            {unknownBranch && <span className="tag alt">not on GitHub</span>}
          </label>
          {branches.length ? (
            <select id="f-GIT_BRANCH" value={branch}
                    onChange={(e) => onChange('GIT_BRANCH', e.target.value)}>
              {unknownBranch && <option value={branch}>{branch} (not on GitHub)</option>}
              {branches.map((b) => <option key={b} value={b}>{b}</option>)}
            </select>
          ) : (
            <input id="f-GIT_BRANCH" value={branch}
                   onChange={(e) => onChange('GIT_BRANCH', e.target.value)} />
          )}
          <code className="pname">GIT_BRANCH</code>
          <p className="help">
            {refs?.repo ? `Branches of ${refs.repo}.` : 'Branch of the sushi repo to run.'}
          </p>
        </div>

        <div className="field">
          <label htmlFor="f-USE_LATEST">Use latest commit on the branch</label>
          <input id="f-USE_LATEST" type="checkbox" checked={useLatest}
                 onChange={(e) => onChange('USE_LATEST', e.target.checked)} />
          <code className="pname">USE_LATEST</code>
          <p className="help">
            {useLatest
              ? head
                ? `Will run ${head.sha} — ${head.subject}`
                : 'Jenkins resolves the branch head at build time.'
              : 'Off: pick the exact commit to run.'}
          </p>
        </div>

        {!useLatest && (
          <div className="field">
            <label htmlFor="f-COMMIT_ID">
              Commit <span className="req">required</span>
            </label>
            {commits.length ? (
              <select id="f-COMMIT_ID" value={values.COMMIT_ID || ''}
                      onChange={(e) => onChange('COMMIT_ID', e.target.value)}>
                <option value="">Select a commit…</option>
                {commits.map((c) => (
                  <option key={c.full_sha} value={c.sha}>
                    {c.sha} · {c.date?.slice(0, 10)} · {c.subject}
                  </option>
                ))}
              </select>
            ) : (
              <input id="f-COMMIT_ID" value={values.COMMIT_ID || ''}
                     onChange={(e) => onChange('COMMIT_ID', e.target.value)} />
            )}
            <code className="pname">COMMIT_ID</code>
            <p className="help">
              {commits.length
                ? `${commits.length} most recent commits on ${branch}.`
                : 'Short SHA of the commit to check out.'}
            </p>
          </div>
        )}
      </div>
    </>
  )
}

function Field({ spec, value, stock, fromPreset, onChange }) {
  const id = `f-${spec.name}`
  const mark = changedClass(value, stock, fromPreset)
  const suggest = useContext(SuggestContext)?.[spec.name]
  const listId = suggest?.values?.length ? `sug-${spec.name}` : undefined
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
      ) : spec.multiline ? (
        // A textarea, not an input, because an input cannot wrap and these
        // values are long enough that a single line hides most of them. Sized
        // from the content so the whole value is visible without dragging;
        // whether a field is multiline comes from the spec, so the element type
        // never changes under a cursor mid-edit.
        <textarea id={id} value={value ?? ''} spellCheck={false}
                  rows={Math.min(6, Math.max(2, Math.ceil(String(value ?? '').length / 46)))}
                  onChange={(e) => onChange(spec.name, e.target.value)} />
      ) : (
        <>
          <input id={id} list={listId} value={value ?? ''} autoComplete="off"
                 onChange={(e) => onChange(spec.name, e.target.value)} />
          {listId && (
            <datalist id={listId}>
              {suggest.values.map((v) => <option key={v} value={v} />)}
            </datalist>
          )}
        </>
      )}
      <code className="pname">{spec.name}</code>
      {spec.help && <p className="help">{spec.help}</p>}
      {suggest?.values?.length > 0 && (
        <p className="help">
          From cellDB: {suggest.values.join(', ')}. Any other value is treated as a CSV file in
          the build directory.
        </p>
      )}
      {suggest?.error && <p className="help warn-text">{suggest.error}</p>}
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
