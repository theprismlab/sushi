import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// In dev, proxy /api to the FastAPI backend so the app is same-origin and
// there is nothing to configure in the client.
export default defineConfig({
  plugins: [react()],
  server: {
    proxy: { '/api': process.env.API_TARGET || 'http://localhost:8000' },
  },
})
