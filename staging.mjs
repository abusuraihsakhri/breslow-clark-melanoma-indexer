const CUTANEOUS_SITES = new Set(["trunk", "extremity", "head_neck", "acral"]);

const SITE_WEIGHTS = {
  trunk: 0.10,
  extremity: 0.00,
  head_neck: 0.30,
  acral: 0.20,
};

function finiteNumber(value, label) {
  const n = Number(value);
  if (!Number.isFinite(n)) throw new Error(`${label} must be a finite number.`);
  return n;
}

export function assignTStage(depth, ulcerated, isInSitu = false) {
  const t = finiteNumber(depth, "Breslow depth");
  if (t < 0) throw new Error("Breslow depth cannot be negative.");
  if (isInSitu && t > 0) {
    throw new Error("In-situ melanoma cannot have a positive Breslow depth.");
  }
  if (isInSitu || t === 0) {
    return { category: "Tis", basis: "Melanoma in situ / Breslow depth 0.0 mm" };
  }
  if (t < 0.8) {
    return {
      category: ulcerated ? "pT1b" : "pT1a",
      basis: `Breslow thickness <0.8 mm; ulceration ${ulcerated ? "present" : "absent"}`,
    };
  }
  if (t <= 1.0) {
    return { category: "pT1b", basis: "Breslow thickness 0.8–1.0 mm" };
  }
  if (t <= 2.0) {
    return {
      category: ulcerated ? "pT2b" : "pT2a",
      basis: `Breslow thickness >1.0–2.0 mm; ulceration ${ulcerated ? "present" : "absent"}`,
    };
  }
  if (t <= 4.0) {
    return {
      category: ulcerated ? "pT3b" : "pT3a",
      basis: `Breslow thickness >2.0–4.0 mm; ulceration ${ulcerated ? "present" : "absent"}`,
    };
  }
  return {
    category: ulcerated ? "pT4b" : "pT4a",
    basis: `Breslow thickness >4.0 mm; ulceration ${ulcerated ? "present" : "absent"}`,
  };
}

export function marginGuidance(depth, isInSitu = false) {
  const t = finiteNumber(depth, "Breslow depth");
  if (isInSitu || t === 0) {
    return { margin: "0.5–1.0 cm", basis: "General clinical margin range for melanoma in situ" };
  }
  if (t <= 1.0) return { margin: "1.0 cm", basis: "General margin for invasive melanoma ≤1.0 mm" };
  if (t <= 2.0) return { margin: "1.0–2.0 cm", basis: "General margin range for melanoma >1.0–2.0 mm" };
  return { margin: "2.0 cm", basis: "General margin for melanoma >2.0 mm" };
}

export function exploratorySlnEstimate({ depth, age, ulcerated, lvi, site }) {
  const t = finiteNumber(depth, "Breslow depth");
  const a = finiteNumber(age, "Age");
  const siteWeight = SITE_WEIGHTS[site] ?? 0;
  const logit = -2.80 + (0.55 * t) + (0.90 * Number(Boolean(ulcerated)))
    + (-0.30 * ((a - 50) / 10)) + (0.80 * Number(Boolean(lvi))) + siteWeight;
  const probability = 1 / (1 + Math.exp(-logit));
  return Math.round(probability * 10000) / 100;
}

export function slnbGuidance({ depth, age, ulcerated, lvi, mitotic, microsatellites }) {
  if (depth === 0) return "Generally not recommended for melanoma in situ";
  if (microsatellites) return "Specialist regional staging review required; do not rely on this prompt alone";
  if (depth < 0.8 && !ulcerated) {
    if (mitotic >= 2 || lvi || age < 40) return "Discuss and consider SLNB in clinical context";
    return "Generally not recommended for routine pT1a cases";
  }
  if ((depth < 0.8 && ulcerated) || (depth >= 0.8 && depth <= 1.0)) {
    return "Discuss and consider SLNB";
  }
  return "Generally recommended when clinically appropriate";
}

export function analyzeCase(input) {
  const specimenId = String(input.specimenId ?? "").trim();
  if (!specimenId) throw new Error("Specimen ID cannot be empty.");

  const age = finiteNumber(input.age, "Age");
  if (age < 0 || age > 120) throw new Error("Age must be between 0 and 120 years.");

  const depth = finiteNumber(input.depth, "Breslow depth");
  const mitotic = finiteNumber(input.mitotic ?? 0, "Mitotic rate");
  if (mitotic < 0) throw new Error("Mitotic rate cannot be negative.");

  const site = String(input.site ?? "");
  if (!CUTANEOUS_SITES.has(site)) throw new Error("Site must be a supported cutaneous melanoma site.");

  const clark = input.clark === "" || input.clark == null ? null : Number(input.clark);
  if (clark !== null && ![1, 2, 3, 4, 5].includes(clark)) throw new Error("Clark level must be I–V.");

  const isInSitu = Boolean(input.isInSitu);
  const ulcerated = Boolean(input.ulcerated);
  if (isInSitu && depth > 0) throw new Error("In-situ melanoma cannot have a positive Breslow depth.");

  const t = assignTStage(depth, ulcerated, isInSitu);
  const margin = marginGuidance(depth, isInSitu);

  const adverse = [];
  if (ulcerated) adverse.push("Ulceration present");
  if (input.lvi) adverse.push("Lymphovascular invasion present");
  if (input.perineural) adverse.push("Perineural / neurotropic invasion present");
  if (input.microsatellites) adverse.push("Microsatellitosis present");
  if (mitotic >= 1) adverse.push(`Mitotic activity ${mitotic.toFixed(1)}/mm²`);
  if (input.regression) adverse.push("Histologic regression reported");

  const alerts = [];
  if (input.microsatellites) {
    alerts.push("Microsatellitosis is a regional metastatic feature; final N category and stage group require complete nodal and distant-metastasis data.");
  }
  if (t.category.startsWith("pT4")) alerts.push("Thick primary melanoma (>4.0 mm): ensure complete guideline-based staging assessment.");
  if (input.perineural) alerts.push("Perineural/neurotropic invasion reported: correlate with subtype and local-management planning.");
  if (clark === 1 && depth > 0) alerts.push("Clark level I is inconsistent with a positive Breslow depth.");

  return {
    specimenId,
    age,
    site,
    depth,
    clark,
    tStage: t,
    margin,
    slnbGuidance: slnbGuidance({
      depth,
      age,
      ulcerated,
      lvi: Boolean(input.lvi),
      mitotic,
      microsatellites: Boolean(input.microsatellites),
    }),
    exploratorySlnProbabilityPct: exploratorySlnEstimate({
      depth,
      age,
      ulcerated,
      lvi: Boolean(input.lvi),
      site,
    }),
    adverse,
    alerts,
  };
}
