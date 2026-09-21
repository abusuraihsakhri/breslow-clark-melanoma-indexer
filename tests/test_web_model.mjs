import assert from "node:assert/strict";
import {
  analyzeCase,
  assignTStage,
  marginGuidance,
} from "../staging.mjs";

assert.equal(assignTStage(0, false, true).category, "Tis");
assert.equal(assignTStage(0.79, false).category, "pT1a");
assert.equal(assignTStage(0.79, true).category, "pT1b");
assert.equal(assignTStage(0.8, false).category, "pT1b");
assert.equal(assignTStage(1.0, true).category, "pT1b");
assert.equal(assignTStage(1.01, false).category, "pT2a");
assert.equal(assignTStage(2.01, true).category, "pT3b");
assert.equal(assignTStage(4.01, false).category, "pT4a");

assert.throws(() => assignTStage(0.8, false, true), /positive Breslow/);
assert.equal(marginGuidance(0.6).margin, "1.0 cm");
assert.equal(marginGuidance(1.5).margin, "1.0–2.0 cm");

const result = analyzeCase({
  specimenId: "WEB-1",
  age: 58,
  site: "trunk",
  depth: 1.45,
  clark: 4,
  mitotic: 3,
  isInSitu: false,
  ulcerated: true,
  lvi: false,
  perineural: false,
  microsatellites: false,
  regression: false,
});
assert.equal(result.tStage.category, "pT2b");
assert.equal(result.clark, 4);
assert.ok(result.exploratorySlnProbabilityPct > 0);

assert.throws(() => analyzeCase({
  specimenId: "WEB-2",
  age: 50,
  site: "mucosal",
  depth: 1,
  mitotic: 0,
}), /cutaneous melanoma site/);

console.log("Browser staging model tests passed.");
