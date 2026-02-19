export default function LoadingSpinner({ message = 'Analyzing genetic data...' }) {
  return (
    <div className="spinner-wrap">
      <div className="spinner"></div>
      <p>{message}</p>
    </div>
  )
}