import { createContext, useContext, useEffect, useState } from 'react'
import { NavLink, Outlet } from 'react-router-dom'
import { api } from './api.js'

const CatalogContext = createContext(null)
export const useCatalog = () => useContext(CatalogContext)

/** Who Jenkins says you are. Identity is not self-reported: the same
 *  credentials that name you are the ones that authorize the build trigger. */
const IdentityContext = createContext([null, () => {}])
export const useIdentity = () => useContext(IdentityContext)

export default function App() {
  const [catalog, setCatalog] = useState(null)
  const [health, setHealth] = useState(null)
  const [identity, setIdentity] = useState(null)
  // Distinct from "signed out": until the cookie check returns we do not know,
  // and showing the wall in the meantime flashes it on every reload.
  const [checked, setChecked] = useState(false)
  const [error, setError] = useState(null)

  useEffect(() => {
    api.session()
      .then((s) => setIdentity(s.identity))
      .catch(() => {})
      .finally(() => setChecked(true))
  }, [])

  // The catalog is only needed once past the wall.
  useEffect(() => {
    if (!identity) return
    api.catalog().then(setCatalog).catch((e) => setError(e.message))
    api.health().then(setHealth).catch(() => {})
  }, [identity])

  if (!checked) return <div className="page"><div className="muted">Loading…</div></div>

  if (!identity) return <SignInWall onSignedIn={setIdentity} />

  if (error) return <div className="page"><div className="banner error">Cannot reach the backend: {error}</div></div>
  if (!catalog) return <div className="page"><div className="muted">Loading…</div></div>

  // Only one of these two directions can break a launch: a parameter we send
  // that the job does not declare makes Jenkins reject the whole build. A
  // parameter the job declares that we never send just keeps its job default.
  const rejects = health?.params_we_send_that_job_rejects || []
  const missing = health?.job_params_not_in_ui || []

  return (
    <CatalogContext.Provider value={catalog}>
      <IdentityContext.Provider value={[identity, setIdentity]}>
        <header className="topbar">
          <div className="brand">
            sushi <span className="muted">PRISM pipeline</span>
            {health?.version && (
              <span className="muted small" title={`branch ${health.version.branch}`}>
                {health.environment} · {health.version.commit}
              </span>
            )}
          </div>
          <nav>
            <NavLink to="/launch">Launch</NavLink>
            <NavLink to="/runs">Run history</NavLink>
            <NavLink to="/screen-types">Screen types</NavLink>
          </nav>
          <div className="who">
            <span className="muted">signed in as</span>
            <strong>{identity.full_name}</strong>
            <button onClick={() => api.signOut().then(() => setIdentity(null))}>Sign out</button>
          </div>
        </header>

        {health && health.environment !== 'production' && (
          <div className="banner warn">
            <strong>{health.environment} instance — not a sandbox</strong>
            There is only one pipeline Jenkins job, so anything you launch here is a real
            pipeline run that writes into the real build directory. Only the run history
            and the screen-type defaults are separate.
          </div>
        )}
        {health && !health.jenkins_reachable && (
          <div className="banner error">
            Jenkins at {health.jenkins_url} is unreachable ({health.jenkins_error}). You can browse
            history but launching will fail.
          </div>
        )}
        {rejects.length > 0 && (
          <div className="banner error">
            The Jenkins job does not declare {rejects.join(', ')}. Launches will be rejected until{' '}
            <code>params.yml</code> and <code>make_config_file.groovy</code> agree.
          </div>
        )}
        {missing.length > 0 && (
          <div className="banner warn">
            The Jenkins job declares parameters this UI does not expose: {missing.join(', ')}.
            Launching still works — those keep whatever default the job defines.
          </div>
        )}

        <main className="page"><Outlet /></main>
      </IdentityContext.Provider>
    </CatalogContext.Provider>
  )
}

function SignInWall({ onSignedIn }) {
  const [user, setUser] = useState('')
  const [password, setPassword] = useState('')
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState(null)

  async function submit(e) {
    e.preventDefault()
    setBusy(true)
    setError(null)
    try {
      onSignedIn(await api.signIn(user.trim(), password.trim()))
    } catch (err) {
      setError(typeof err.detail === 'string' ? err.detail : err.message)
      setBusy(false)
    }
  }

  return (
    <div className="wall">
      <form className="card wall-card" onSubmit={submit}>
        <div className="brand">sushi <span className="muted">PRISM pipeline</span></div>
        <p className="muted">
          Sign in with your Jenkins credentials. Jenkins requires them to queue a build, and they
          are what the run log records as the launcher.
        </p>

        <div className="field">
          <label htmlFor="wall-user">Jenkins user ID</label>
          <input id="wall-user" value={user} onChange={(e) => setUser(e.target.value)}
                 autoComplete="username" autoFocus />
        </div>
        <div className="field">
          <label htmlFor="wall-password">Password</label>
          <input id="wall-password" type="password" value={password}
                 onChange={(e) => setPassword(e.target.value)} autoComplete="current-password" />
          <p className="help">
            If Jenkins rejects it, an API token from <code>/user/&lt;you&gt;/configure</code> works
            here too — some Jenkins configurations only accept tokens over the API.
          </p>
        </div>

        {error && <div className="banner error">{error}</div>}

        <button className="primary" disabled={busy || !user.trim() || !password.trim()}>
          {busy ? 'Checking…' : 'Sign in'}
        </button>
      </form>
    </div>
  )
}
