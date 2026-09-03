import { createContext, useContext, useEffect, useState } from 'react'
import { NavLink, Outlet } from 'react-router-dom'
import { api } from './api.js'

const CatalogContext = createContext(null)
export const useCatalog = () => useContext(CatalogContext)

/** There is no auth, so identity is self-reported and remembered locally.
 *  It is still the only record of who launched a build, so we require it. */
export function useIdentity() {
  const [who, setWho] = useState(() => localStorage.getItem('sushi.who') || '')
  useEffect(() => { localStorage.setItem('sushi.who', who) }, [who])
  return [who, setWho]
}

export default function App() {
  const [catalog, setCatalog] = useState(null)
  const [health, setHealth] = useState(null)
  const [error, setError] = useState(null)
  const [who, setWho] = useIdentity()

  useEffect(() => {
    api.catalog().then(setCatalog).catch((e) => setError(e.message))
    api.health().then(setHealth).catch(() => {})
  }, [])

  if (error) return <div className="page"><div className="banner error">Cannot reach the backend: {error}</div></div>
  if (!catalog) return <div className="page"><div className="muted">Loading…</div></div>

  const drift = [
    ...(health?.params_we_send_that_job_rejects || []),
    ...(health?.job_params_not_in_ui || []),
  ]

  return (
    <CatalogContext.Provider value={catalog}>
      <header className="topbar">
        <div className="brand">
          sushi <span className="muted">PRISM pipeline</span>
        </div>
        <nav>
          <NavLink to="/launch">Launch</NavLink>
          <NavLink to="/runs">Run history</NavLink>
        </nav>
        <label className="who">
          <span className="muted">you are</span>
          <input
            value={who}
            onChange={(e) => setWho(e.target.value)}
            placeholder="your name"
            aria-label="Your name, recorded against every run you launch"
          />
        </label>
      </header>

      {health && !health.jenkins_reachable && (
        <div className="banner error">
          Jenkins at {health.jenkins_url} is unreachable ({health.jenkins_error}). You can browse
          history but launching will fail.
        </div>
      )}
      {drift.length > 0 && (
        <div className="banner warn">
          Parameter drift between this UI and the Jenkins job: {drift.join(', ')}. Launches may be
          rejected until <code>params.yml</code> and <code>make_config_file.groovy</code> agree.
        </div>
      )}

      <main className="page"><Outlet /></main>
    </CatalogContext.Provider>
  )
}
