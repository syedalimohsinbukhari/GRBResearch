# Model Selection Criteria

This document outlines the criteria used for selecting the best GRB spectral models from SAFE to GOOD to BEST classifications.

---

## Table of Contents

1. [Model Complexity Hierarchy](#model-complexity-hierarchy)
2. [C-Stat Comparison Logic](#c-stat-comparison-logic)
3. [Error Criteria](#error-criteria)
4. [Model Classification Flow](#model-classification-flow)
5. [Model Hierarchy Analysis](#model-hierarchy-analysis)

---

## Model Complexity Hierarchy

### Single (Base) Models

Ordered by complexity:

| Model | Free Params | Complexity Order |
|-------|-------------|------------------|
| PL    | 2           | 0                |
| CPL   | 3           | 1                |
| BAND  | 4           | 2                |
| SBPL  | 4           | 2                |

### Additional Components

Each component adds free parameters:

| Component | Additional Free Params | Description          |
|-----------|------------------------|----------------------|
| +BB       | +2                     | Blackbody component  |
| +PL       | +2                     | Power-law component  |

### Model Groups

Models are organized into four groups for selection:

1. **BASE**: `PL`, `CPL`, `BAND`, `SBPL`
2. **+BB**: `PL_BB`, `CPL_BB`, `BAND_BB`, `SBPL_BB`
3. **+PL**: `CPL_PL`, `BAND_PL`, `SBPL_PL`
4. **+PL+BB**: `CPL_PL_BB`, `BAND_PL_BB`, `SBPL_PL_BB`

---

## C-Stat Comparison Logic

The `compare_models()` function determines which model is better based on C-Stat values and model complexity.

### 1. NaN Handling

- If **one model** has NaN c-stat → choose the **valid model**
- If **both models** have NaN → choose the **simpler model** (by base model order)

### 2. Equal Free Parameters

When both models have the same number of free parameters:

- If c-stat values are **equal** → choose the **simpler model** (by `complexity_key`)
- Otherwise → choose the model with **lower c-stat**

### 3. Different Free Parameters (Key Rule)

When models have different numbers of free parameters:

```text
required_improvement = 9 × (complex_free_params - simple_free_params)
actual_improvement   = simple_cstat - complex_cstat
```

**Decision Rule:**
- Choose **complex model** if `actual_improvement >= required_improvement`
- Otherwise → choose **simpler model**

**Example:**
- Simple model: PL (2 params), c-stat = 500
- Complex model: BAND (4 params), c-stat = 482
- Required improvement: 9 × (4 - 2) = 18
- Actual improvement: 500 - 482 = 18
- **Result**: Choose BAND (improvement meets threshold)

---

## Error Criteria

Models must pass error criteria to be considered **SAFE**. The `model_passes_error_criteria()` function evaluates each parameter.

### Failure Conditions

A model **FAILS** if any parameter satisfies:

| Condition | Result | Reason |
|-----------|--------|--------|
| `value == 0` **AND** `error != 0` | **FAIL** | Parameter is zero but has non-zero error |
| `abs(error) >= limit` | **FAIL** | Error is too large relative to value |

### Error Limits by Parameter Type

| Parameter Type | Condition | Limit Formula | Notes |
|----------------|-----------|---------------|-------|
| **Default** | All other params | `0.4 × abs(value)` | 40% relative error |
| **index2** (BAND/SBPL) | `loose_criteria=True` | `1.0 × abs(value)` | 100% relative error |
| **index2** (BAND/SBPL) | `loose_criteria=False` | `0.7 × abs(value)` | 70% relative error |
| **index2_pl** (+PL models) | `loose_criteria=True` | `1.0 × abs(value)` | 100% relative error |
| **index2_pl** (+PL models) | `loose_criteria=False` | `0.7 × abs(value)` | 70% relative error |

### Parameters

- **`par_constraint`**: Default = `0.4` (40% relative error allowed)
- **`loose_criteria`**: Default = `False` (stricter error limits on `index2` parameters)

---

## Model Classification Flow

### 1. SAFE Models

**Definition**: Models that pass error criteria

**Selection Process:**
1. Check if model `.fit` file exists
2. Read parameter values and errors
3. Apply error criteria to each parameter
4. Model is **SAFE** if all parameters pass

**Function**: `filter_models_by_error(c_stats, folder_path, candidates, **kwargs)` (`safe_good_best.py`) — renamed from `list_safe_models`, which no longer exists under that name.

### 2. GOOD Models

**Definition**: Best model per group among SAFE models

**Selection Process:**

For each group (`BASE`, `BB`, `PL`, `PLBB`):

1. **Filter**: Keep only SAFE models in the group
2. **Compare**: Use `compare_models()` to find the best within the group
3. **Result**: One GOOD model per group (if any exist)

**Groups:**
- `BASE`: Best among `PL`, `CPL`, `BAND`, `SBPL`
- `BB`: Best among `PL_BB`, `CPL_BB`, `BAND_BB`, `SBPL_BB`
- `PL`: Best among `CPL_PL`, `BAND_PL`, `SBPL_PL`
- `PLBB`: Best among `CPL_PL_BB`, `BAND_PL_BB`, `SBPL_PL_BB`

**Function**: `pick_best_model(c_stats, candidates, group_name, folder_path=None, is_separate_group=-1, **kwargs)` (`safe_good_best.py`), which internally calls `filter_models_by_error` then `pick_best_in_group` — renamed from `compute_good_models`, which no longer exists under that name.

### 3. BEST Model

**Definition**: Final selected model from model hierarchy analysis

**Selection Process:**
1. Analyze GOOD models using `analyze_model_hierarchy()`
2. Apply comparison rules based on component additions
3. Assign status: `ACCEPTED`, `REJECTED`, `UNNECESSARY`, or `INVALID`
4. The **BEST** model is the one marked as `ACCEPTED`

---

## Model Hierarchy Analysis

The `analyze_model_hierarchy()` function determines which GOOD model should be the BEST.

### Comparison Rules

The comparison is always relative to the **BASE** model:

#### Rule 1: Single Component Addition (BASE → BASE_XX)

**Threshold**: `Δc-stat > 25`

- Example: `CPL` → `CPL_BB`
- If `c-stat(CPL) - c-stat(CPL_BB) > 25` → **ACCEPT** `CPL_BB`
- Otherwise → **REJECT** `CPL_BB`

#### Rule 2: Double Component Addition (BASE → BASE_XX_YY)

Two scenarios depending on single component results:

##### Scenario A: At least one single extension is ACCEPTED

**Threshold**: `Δc-stat > 25` (from best single extension)

- Example: `CPL_BB` is ACCEPTED, evaluating `CPL_PL_BB`
- Compare against best of `CPL_BB` or `CPL_PL`
- If `c-stat(best_single) - c-stat(CPL_PL_BB) > 25` → **ACCEPT**
- Otherwise → **REJECT**

##### Scenario B: Both single extensions are REJECTED

**Threshold**: `Δc-stat > 50` (from BASE)

- Example: Both `CPL_BB` and `CPL_PL` are REJECTED, evaluating `CPL_PL_BB`
- If `c-stat(CPL) - c-stat(CPL_PL_BB) > 50` → **ACCEPT**
- Otherwise → **REJECT**

### Model Status Flags

| Status | Meaning |
|--------|---------|
| **ACCEPTED** | Model is significantly better and should be used |
| **REJECTED** | Model does not provide sufficient improvement |
| **UNNECESSARY** | Model is beaten by a more complex variant |
| **INVALID** | Model does not contain the BASE model (inconsistent comparison) |

### BASE Model Status

The BASE model receives special handling:

- If **any** model beats BASE → BASE is marked **UNNECESSARY**
- If **no** model beats BASE → BASE is marked **ACCEPTED**

### Summary Table

| BASE Extensions | Improvement Threshold | Comparison Reference |
|-----------------|----------------------|----------------------|
| 1 component (+BB or +PL) | > 25 | BASE model |
| 2 components (+PL+BB) with 1 single ACCEPTED | > 25 | Best single extension |
| 2 components (+PL+BB) with both singles REJECTED | > 50 | BASE model |

**Fallback for models that don't decompose into a clean single/double extension** (a model name containing the BASE name but not matching either category — `analyze_model_hierarchy`'s final loop, `grb_utils.py`): uses `cstat_threshold(delta_k=2*additional_components, sigma=5.0)` instead of the flat 25/50 above — the same 5σ, χ²(Δk) threshold used elsewhere in the paper (28.74 for Δk=2, 36.86 for Δk=4), not yet reconciled with this flat-25/50 rule for the primary categorized path.

---

## Example Workflow

### Scenario: Evaluating models for an epoch

**Available models**: `CPL`, `CPL_BB`, `CPL_PL`, `CPL_PL_BB`

**C-Stat values**:
- `CPL`: 500
- `CPL_BB`: 460
- `CPL_PL`: 490
- `CPL_PL_BB`: 430

#### Step 1: SAFE Classification
All models pass error criteria → All are **SAFE**

#### Step 2: GOOD Classification
- `BASE` group: `CPL` is the only candidate → `CPL` is **GOOD**
- `BB` group: `CPL_BB` is the only candidate → `CPL_BB` is **GOOD**
- `PL` group: `CPL_PL` is the only candidate → `CPL_PL` is **GOOD**
- `PLBB` group: `CPL_PL_BB` is the only candidate → `CPL_PL_BB` is **GOOD**

#### Step 3: Hierarchy Analysis

1. **Evaluate** `CPL_BB`:
   - Improvement: 500 - 460 = 40 > 25 ✓
   - Status: **ACCEPTED**

2. **Evaluate** `CPL_PL`:
   - Improvement: 500 - 490 = 10 < 25 ✗
   - Status: **REJECTED**

3. **Evaluate** `CPL_PL_BB`:
   - One single extension (`CPL_BB`) is ACCEPTED
   - Compare against `CPL_BB`: 460 - 430 = 30 > 25 ✓
   - Status: **ACCEPTED**

4. **Evaluate BASE** (`CPL`):
   - `CPL_BB` and `CPL_PL_BB` are ACCEPTED
   - Status: **UNNECESSARY**

#### Final Result:
**BEST** model: `CPL_PL_BB` (most complex ACCEPTED model)

---

## Notes

- All comparisons prioritize **simpler models** unless complexity provides statistically significant improvement
- The threshold of **9 c-stat units per additional free parameter** ensures proper penalization of complexity
- Error criteria prevent over-fitted models with poorly constrained parameters from being selected
- The hierarchy analysis ensures that component additions are justified by sufficient improvement in fit quality
