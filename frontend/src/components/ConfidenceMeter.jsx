export default function ConfidenceMeter({ score, riskLabel }) {
  const percentage = Math.round((score || 0) * 100)

  const getColor = () => {
    if (riskLabel === 'Safe') return '#8C7B58'
    if (riskLabel === 'Adjust Dosage') return '#DAA99A'
    if (riskLabel === 'Toxic') return '#C75C5F'
    return '#B7AFA7'
  }

  const color = getColor()

  // Semi-circle math
  const radius = 54
  const circumference = Math.PI * radius
  const offset = circumference - (percentage / 100) * circumference

  return (
    <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center', gap: 8 }}>
      <div className="micro" style={{ fontSize: 10, letterSpacing: 2, fontFamily: 'var(--mono)', color: 'var(--text-dim)', textTransform: 'uppercase' }}>
        Clinical Confidence
      </div>

      <div style={{ position: 'relative', width: 120, height: 64, overflow: 'hidden' }}>
        <svg width="120" height="120" style={{ position: 'absolute', top: 0, left: 0 }}>
          {/* Background track */}
          <circle
            cx="60" cy="60" r={radius}
            fill="none"
            stroke="var(--border2)"
            strokeWidth="10"
            strokeDasharray={circumference}
            strokeDashoffset={0}
            strokeLinecap="round"
            transform="rotate(-180 60 60)"
            style={{ opacity: 0.4 }}
          />
          {/* Foreground arc */}
          <circle
            cx="60" cy="60" r={radius}
            fill="none"
            stroke={color}
            strokeWidth="10"
            strokeDasharray={circumference}
            strokeDashoffset={offset}
            strokeLinecap="round"
            transform="rotate(-180 60 60)"
            style={{ transition: 'stroke-dashoffset 1s ease-in-out, stroke 0.5s ease' }}
          />
        </svg>

        {/* Center label */}
        <div style={{
          position: 'absolute',
          bottom: 0,
          width: '100%',
          textAlign: 'center',
        }}>
          <div style={{ fontFamily: 'var(--mono)', fontSize: 20, fontWeight: 700, color: color }}>
            {percentage}%
          </div>
        </div>
      </div>

      {/* Risk label pill */}
      <span style={{
        padding: '3px 12px',
        borderRadius: 100,
        fontSize: 10,
        fontFamily: 'var(--mono)',
        fontWeight: 700,
        letterSpacing: 1,
        color: color,
        border: `1.5px solid ${color}`,
        background: `${color}22`,
      }}>
        {riskLabel || 'UNKNOWN'}
      </span>
    </div>
  )
}