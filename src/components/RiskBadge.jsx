const RISK_CLASS = {
  'Safe': 'Safe',
  'Adjust Dosage': 'Adjust',
  'Toxic': 'Toxic',
  'Ineffective': 'Ineffective',
  'Unknown': 'Unknown',
}

export default function RiskBadge({ riskLabel }) {
  const cls = RISK_CLASS[riskLabel] || 'Unknown'
  return (
    <span className={'risk-badge ' + cls}>
      {riskLabel || 'Unknown'}
    </span>
  )
}