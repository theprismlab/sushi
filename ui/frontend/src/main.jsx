import React from 'react'
import ReactDOM from 'react-dom/client'
import { createBrowserRouter, RouterProvider, Navigate } from 'react-router-dom'
import App from './App.jsx'
import Launch from './Launch.jsx'
import History from './History.jsx'
import RunDetail from './RunDetail.jsx'
import ScreenTypes from './ScreenTypes.jsx'
import './styles.css'

const router = createBrowserRouter([
  {
    path: '/',
    element: <App />,
    children: [
      { index: true, element: <Navigate to="/launch" replace /> },
      { path: 'launch', element: <Launch /> },
      { path: 'runs', element: <History /> },
      { path: 'runs/:id', element: <RunDetail /> },
      { path: 'screen-types', element: <ScreenTypes /> },
    ],
  },
])

ReactDOM.createRoot(document.getElementById('root')).render(
  <React.StrictMode>
    <RouterProvider router={router} />
  </React.StrictMode>,
)
