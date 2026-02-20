# 🧬 PharmaGuard — Clinical AI for Pharmacogenomic Risk Intelligence

**PharmaGuard** is an end-to-end pharmacogenomics decision-support platform built during a high-pressure hackathon to transform raw genomic data into fast, explainable, clinician-ready insights.

While many teams focused only on AI generation or UI dashboards, **we engineered a clinically grounded system** — combining genomics parsing, CPIC-aligned risk modeling, and safety-focused AI explanations into a single real-time workflow.

---

## 🚀 What PharmaGuard Does

Upload a patient **VCF file** + select a drug → PharmaGuard instantly delivers:

* 🧬 Gene + Diplotype interpretation
* ⚠️ CPIC-aligned risk classification
* 📊 Pharmacogenomic profile with detected variants
* 🤖 Grounded AI clinical explanation (mechanism-aware)
* 🔥 Risk Heatmap visualization
* 💬 “Ask PharmaGuard” clinician assistant

All powered by a **privacy-first local LLM pipeline**.

---

## 🏆 Why PharmaGuard Stands Out (Even Among Similar PS Builds)

Many hackathon solutions stopped at *“AI generates medical text.”*
PharmaGuard goes further — focusing on **clinical reliability, explainability, and real workflow design**.

### ✅ 1. Clinically Grounded AI — Not Generic LLM Output

* AI responses are constrained by structured pharmacogenomic data
* Variant citations and biological mechanisms included
* CPIC-aligned reasoning enforced
* Hallucination-resistant prompt architecture

### ⚡ 2. Real-Time Non-Blocking Architecture

* Risk engine responds instantly
* LLM explanation loads asynchronously
* Ultra-Turbo async pipeline reduces perceived latency
* Model warm-keep and token-limited inference for speed

### 🔥 3. Risk Heatmap Mode (Hackathon Differentiator)

Dynamic UI highlighting toxicity intensity using genomic risk signals — translating complex PGx data into an immediately interpretable clinical visual.

### 💬 4. Ask PharmaGuard — Context-Aware Clinical AI

Unlike basic chatbots:

* Grounded in current patient genotype + drug
* Mechanism-focused answers
* Clinician-tone generation
* Safety disclaimers enforced

### 🛡️ 5. Safety-First Design Philosophy

* Structured schema validation
* Fallback explanation logic
* No hallucinated genes or variants allowed
* Built for **decision support**, not autonomous prescribing

---

## 🧠 Tech Stack

**Backend**

* FastAPI
* Pydantic
* Ollama (Llama3 local model)
* Async pipeline architecture

**Frontend**

* React + Vite
* Real-time state updates
* Clinical heatmap UI system

**AI Layer**

* Local LLM inference
* Structured prompt grounding
* Non-blocking async explanation pipeline

---

## ⚙️ Architecture Highlights

```
VCF Upload → Variant Parser → Risk Engine → Clinical Recommendation
                              ↓
                        Async LLM Explanation
                              ↓
                    Heatmap + Clinician UI
```

---

## 🧪 Key Engineering Achievements During Hackathon

* Built full pharmacogenomic pipeline from scratch
* Integrated local LLM safely into clinical workflow
* Designed hallucination-resistant explanation system
* Implemented real-time UX despite heavy AI inference
* Delivered complete frontend + backend integration under time constraints

---

## 🎥 Demo

👉 Full walkthrough video included with this submission.

---

## 🙌 Acknowledgements

Massive thanks to **RIFT by PWIOI** for hosting an incredible hackathon environment that encouraged deep technical execution, clinical responsibility, and bold experimentation.

---

## ⚠️ Disclaimer

PharmaGuard provides AI-assisted pharmacogenomic insights for educational and clinical support purposes only. Final prescribing decisions must be made by licensed healthcare professionals.

---

## 💥 Vision

PharmaGuard is not just a hackathon project — it’s a step toward making pharmacogenomics **interpretable, actionable, and safe** in real clinical workflows.

Because precision medicine shouldn’t feel like decoding a genome alone.
