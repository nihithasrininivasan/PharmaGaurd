import { useEffect, useState } from 'react'

const GENERATING = 'Generating Clinical Explanation...'

export default function LoadingSpinner({ message = 'Analyzing genetic data...' }) {
  const [dots, setDots] = useState('')

  useEffect(() => {
    if (message !== GENERATING) return
    const interval = setInterval(() => {
      setDots(d => d.length >= 3 ? '' : d + '.')
    }, 400)
    return () => clearInterval(interval)
  }, [message])

  const isGenerating = message === GENERATING

  return (
    <div className="spinner-wrap">
      <div className="spinner"></div>
      {isGenerating ? (
        <p>
          Generating Clinical Explanation
          <span style={{ display: 'inline-block', width: 24, textAlign: 'left' }}>{dots}</span>
        </p>
      ) : (
        <p>{message}</p>
      )}
    </div>
  )
}