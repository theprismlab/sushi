import { useEffect, useMemo, useState } from 'react'
import { api } from './api.js'
import { useCatalog } from './App.jsx'

/** Editor for the screen-type defaults.
 *
 *  These used to be derived by hand from the most recent build of each type.
 *  They are now a config file this page owns: presets.yml in git seeds it, and
 *  the first save writes a store that takes over from then on. */
export default function ScreenTypes() {
  const catalog = useCatalog()
  const [state, setState] = useState(null)
  const [draft, setDraft] = useState(null)
  const [problems, setProblems] = useState([])
  const [busy, setBusy] = useState(false)
  const [saved, setSaved] = useState(false)

  const specs = useMemo(
    () => [...catalog.params].sort((a, b) => a.name.localeCompare(b.name)),
    [catalog],
  )
  const index = useMemo(
    () => Object.fromEntries(catalog.params.map((p) => [p.name, p])),
    [catalog],
  )

  useEffect(() => { load() }, [])

  function load() {
    api.screenTypes().then((d) => { setState(d); setDraft(structuredClone(d.presets)) })
  }

  if (!draft) return <div className="muted">Loading screen types…</div>

  const dirty = JSON.stringify(draft) !== JSON.stringify(state.presets)

  function update(i, patch) {
    setDraft((d) => d.map((p, j) => (j === i ? { ...p, ...patch } : p)))
    setSaved(false)
  }

  function setParam(i, name, value) {
    setDraft((d) => d.map((p, j) => (j === i ? { ...p, params: { ...p.params, [name]: value } } : p)))
    setSaved(false)
  }

  function removeParam(i, name) {
    setDraft((d) => d.map((p, j) => {
      if (j !== i) return p
      const { [name]: _drop, ...rest } = p.params
      return { ...p, params: rest }
    }))
    setSaved(false)
  }

  function addParam(i, name) {
    if (!name) return
    const spec = index[name]
    // Seed with the catalog default so the row starts from a legal value.
    setParam(i, name, spec.type === 'bool' ? !!spec.default : String(spec.default ?? ''))
  }

  function addScreenType() {
    setDraft((d) => [...d, { id: '', label: '', dir: '', reference: '', params: {} }])
    setSaved(false)
  }

  async function save() {
    setBusy(true)
    setProblems([])
    try {
      const d = await api.saveScreenTypes(draft)
      setState(d)
      setDraft(structuredClone(d.presets))
      setSaved(true)
    } catch (e) {
      setProblems(e.detail?.problems || [e.message])
    } finally {
      setBusy(false)
    }
  }

  async function reset() {
    setBusy(true)
    try {
      const d = await api.resetScreenTypes()
      setState(d)
      setDraft(structuredClone(d.presets))
      setProblems([])
    } finally {
      setBusy(false)
    }
  }

  return (
    <div className="launch">
      <section className="card">
        <h2>Screen type defaults</h2>
        <p className="muted">
          What each screen type prefills on the launch form. Only values that differ from the
          catalog default belong here — everything else falls through to <code>params.yml</code>,
          which keeps the diff against stock readable.
        </p>
        <p className="help">
          {state.source === 'store' ? (
            <>Reading from <code>{state.store}</code>. Edits here are live for everyone.</>
          ) : (
            <>
              Reading the shipped <code>{state.shipped}</code>. Saving writes{' '}
              <code>{state.store}</code>, which takes over from then on.
            </>
          )}
        </p>
      </section>

      {draft.map((preset, i) => (
        <section className="card" key={i}>
          <div className="grid">
            <div className="field">
              <label>Id <span className="req">required</span></label>
              <input value={preset.id} onChange={(e) => update(i, { id: e.target.value })}
                     placeholder="MTS_SEQ" />
              <p className="help">Shown on the launch form. Letters, digits, _ and - only.</p>
            </div>
            <div className="field">
              <label>Label <span className="req">required</span></label>
              <input value={preset.label} onChange={(e) => update(i, { label: e.target.value })}
                     placeholder="MTS — Multiplexed Titration Screen" />
            </div>
            <div className="field">
              <label>Directory <span className="req">required</span></label>
              <input value={preset.dir} onChange={(e) => update(i, { dir: e.target.value })}
                     placeholder="MTS_SEQ" />
              <p className="help">
                Under the prismSeq root. Used to infer the screen type from a build path, so it
                should match where builds of this type actually live.
              </p>
            </div>
            <div className="field">
              <label>Reference build</label>
              <input value={preset.reference || ''}
                     onChange={(e) => update(i, { reference: e.target.value })}
                     placeholder="MTS032" />
              <p className="help">Provenance only — where these values were taken from.</p>
            </div>
          </div>

          <h3>Parameter overrides</h3>
          {Object.keys(preset.params || {}).length === 0 ? (
            <p className="help">No overrides: this screen type uses the catalog defaults verbatim.</p>
          ) : (
            <div className="table-scroll">
              <table className="mini full">
                <thead><tr><th>Parameter</th><th>Value</th><th>Catalog default</th><th /></tr></thead>
                <tbody>
                  {Object.entries(preset.params).sort().map(([name, value]) => {
                    const spec = index[name]
                    return (
                      <tr key={name}>
                        <td><code className="wrap">{name}</code></td>
                        <td>
                          {!spec ? (
                            <span className="warn-text">unknown parameter</span>
                          ) : spec.type === 'bool' ? (
                            <input type="checkbox" checked={!!value}
                                   onChange={(e) => setParam(i, name, e.target.checked)} />
                          ) : spec.type === 'choice' ? (
                            <select value={value} onChange={(e) => setParam(i, name, e.target.value)}>
                              {spec.choices.map((c) => <option key={c} value={c}>{c}</option>)}
                            </select>
                          ) : (
                            <input value={value}
                                   onChange={(e) => setParam(i, name, e.target.value)} />
                          )}
                        </td>
                        <td className="muted wrap small">{String(spec?.default ?? '')}</td>
                        <td>
                          <button onClick={() => removeParam(i, name)} title="Remove override">
                            ×
                          </button>
                        </td>
                      </tr>
                    )
                  })}
                </tbody>
              </table>
            </div>
          )}

          <div className="launch-row">
            <select value="" onChange={(e) => addParam(i, e.target.value)} style={{ width: 320 }}>
              <option value="">Add an override…</option>
              {specs
                .filter((s) => !(s.name in (preset.params || {})))
                .map((s) => (
                  <option key={s.name} value={s.name}>{s.name} — {s.label || s.type}</option>
                ))}
            </select>
            <button onClick={() => { setDraft((d) => d.filter((_, j) => j !== i)); setSaved(false) }}>
              Delete {preset.id || 'this screen type'}
            </button>
          </div>
        </section>
      ))}

      <section className="card sticky">
        {problems.length > 0 && (
          <div className="banner error">
            <strong>Not saved</strong>
            <ul>{problems.map((p) => <li key={p}>{p}</li>)}</ul>
          </div>
        )}
        <div className="launch-row">
          <button className="primary" disabled={busy || !dirty} onClick={save}>
            {busy ? 'Saving…' : 'Save screen types'}
          </button>
          <button onClick={addScreenType}>Add a screen type</button>
          <button disabled={busy || state.source !== 'store'} onClick={reset}
                  title="Delete the local store and use the presets.yml from git">
            Revert to shipped defaults
          </button>
          {saved && <span className="muted">Saved.</span>}
          {dirty && !saved && <span className="muted">Unsaved changes.</span>}
        </div>
      </section>
    </div>
  )
}
