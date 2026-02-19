const SUPPORTED_DRUGS = [
  'CODEINE',
  'WARFARIN',
  'CLOPIDOGREL',
  'SIMVASTATIN',
  'AZATHIOPRINE',
  'FLUOROURACIL',
]

export default function DrugInput({ value, onChange }) {
  return (
    <div className="drug-section">
      <p className="section-label">Drug Selection</p>
      <div className="drug-chips">
        {SUPPORTED_DRUGS.map((drug) => (
          <button
            key={drug}
            onClick={() => onChange(drug)}
            className={'drug-chip' + (value === drug ? ' selected' : '')}
          >
            {drug}
          </button>
        ))}
      </div>
      <input
        type="text"
        value={value}
        onChange={(e) => onChange(e.target.value.toUpperCase())}
        placeholder="Or type a drug name..."
        className="drug-input"
      />
    </div>
  )
}