import { analyzeCase } from "./staging.mjs";

const $ = (id) => document.getElementById(id);
const form = $("caseForm");
const workspace = $("workspace");
let lastResult = null;

function setTheme(theme) {
  document.documentElement.dataset.theme = theme;
  localStorage.setItem("theme", theme);
  $("themeToggle").setAttribute("aria-label", theme === "dark" ? "Switch to light mode" : "Switch to dark mode");
}

const savedTheme = localStorage.getItem("theme");
setTheme(savedTheme || (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light"));

$("themeToggle").addEventListener("click", () => {
  setTheme(document.documentElement.dataset.theme === "dark" ? "light" : "dark");
});

function setMobileView(view) {
  workspace.dataset.mobileView = view;
  document.querySelectorAll(".tab-button").forEach((button) => {
    button.classList.toggle("active", button.dataset.view === view);
  });
}

document.querySelectorAll(".tab-button").forEach((button) => {
  button.addEventListener("click", () => setMobileView(button.dataset.view));
});

function syncInSitu() {
  const checked = $("inSitu").checked;
  $("depth").disabled = checked;
  if (checked) $("depth").value = "0";
}

$("inSitu").addEventListener("change", syncInSitu);

function readInput() {
  return {
    specimenId: $("specimenId").value,
    age: $("age").value,
    site: $("site").value,
    depth: $("depth").value,
    clark: $("clark").value,
    mitotic: $("mitotic").value,
    isInSitu: $("inSitu").checked,
    ulcerated: $("ulcerated").checked,
    lvi: $("lvi").checked,
    perineural: $("perineural").checked,
    microsatellites: $("microsatellites").checked,
    regression: $("regression").checked,
  };
}

function fillList(element, items, emptyText) {
  element.replaceChildren();
  const values = items.length ? items : [emptyText];
  for (const item of values) {
    const li = document.createElement("li");
    li.textContent = item;
    element.appendChild(li);
  }
}

function clarkLabel(value) {
  if (!value) return "Not specified";
  return ["", "I", "II", "III", "IV", "V"][value];
}

function renderResult(result) {
  $("emptyState").hidden = true;
  $("resultContent").hidden = false;
  $("tStage").textContent = result.tStage.category;
  $("tBasis").textContent = result.tStage.basis;
  $("clarkResult").textContent = result.clark ? `Level ${clarkLabel(result.clark)}` : "Not specified";
  $("marginResult").textContent = result.margin.margin;
  $("marginBasis").textContent = result.margin.basis;
  $("slnbResult").textContent = result.slnbGuidance;
  $("riskResult").textContent = `${result.exploratorySlnProbabilityPct.toFixed(2)}%`;
  fillList($("adverseList"), result.adverse, "No listed adverse features from the entered fields.");
  fillList($("alertList"), result.alerts, "No additional staging notes.");
  $("downloadBtn").disabled = false;
}

form.addEventListener("submit", (event) => {
  event.preventDefault();
  $("formError").textContent = "";
  try {
    lastResult = analyzeCase(readInput());
    renderResult(lastResult);
    if (matchMedia("(max-width: 900px)").matches) setMobileView("results");
  } catch (error) {
    lastResult = null;
    $("formError").textContent = error instanceof Error ? error.message : "Unable to analyse this case.";
  }
});

$("sampleBtn").addEventListener("click", () => {
  $("specimenId").value = "MEL-DEMO-01";
  $("age").value = "62";
  $("site").value = "trunk";
  $("depth").value = "1.45";
  $("clark").value = "4";
  $("mitotic").value = "3";
  $("inSitu").checked = false;
  $("ulcerated").checked = true;
  $("lvi").checked = false;
  $("perineural").checked = false;
  $("microsatellites").checked = false;
  $("regression").checked = false;
  syncInSitu();
  $("formError").textContent = "";
});

$("resetBtn").addEventListener("click", () => {
  form.reset();
  $("specimenId").value = "MEL-001";
  $("age").value = "55";
  $("site").value = "trunk";
  $("depth").value = "0.75";
  $("mitotic").value = "0";
  syncInSitu();
  lastResult = null;
  $("formError").textContent = "";
  $("emptyState").hidden = false;
  $("resultContent").hidden = true;
  $("downloadBtn").disabled = true;
  setMobileView("input");
});

$("downloadBtn").addEventListener("click", () => {
  if (!lastResult) return;
  const blob = new Blob([JSON.stringify(lastResult, null, 2)], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = `${lastResult.specimenId.replace(/[^a-z0-9_-]+/gi, "_")}_melanoma_index.json`;
  document.body.appendChild(anchor);
  anchor.click();
  anchor.remove();
  URL.revokeObjectURL(url);
});

syncInSitu();
