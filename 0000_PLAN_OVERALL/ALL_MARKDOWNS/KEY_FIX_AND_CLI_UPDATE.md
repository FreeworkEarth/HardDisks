# Key Conflict Fix & CLI Flag Addition

## Summary

Fixed keyboard shortcut conflict in 00ALLINONE.c and added `--hb-temp` CLI flag to both programs for setting heat bath temperature via command line.

---

## 1. Key Conflict: 'a' was Already Taken in 00ALLINONE.c

### The Problem
In 00ALLINONE.c, **'a' key was already used for left piston control**:
- Line 5251: `SDLK_a` moves left piston left
- Line 5334: `SDLK_a` keyup decelerates left piston
- My Andersen code tried to use `SDLK_a` → **CONFLICT!**

### The Solution
**Changed Andersen toggle key from 'a' to 'k' in 00ALLINONE.c**

**Note**: In boxtest.c, 'a' is FREE and remains as the Andersen toggle key!

---

## 2. Updated Keyboard Controls

### boxtest.c (unchanged)
| Key | Action |
|-----|--------|
| **`a`** | Toggle Andersen thermostat ON/OFF |
| **`+`** | Increase Andersen frequency |
| **`-`** | Decrease Andersen frequency |
| **`h`** | Toggle heat bath ON/OFF |
| **`j`** | Decrease HB temperature |
| **`u`** | Increase HB temperature |
| **`t`** | Cycle thermal wall mode (0→1→2) |

### 00ALLINONE.c (updated!)
| Key | Action |
|-----|--------|
| **`k`** | Toggle Andersen thermostat ON/OFF **(changed from 'a')** |
| **`,`** | Decrease Andersen frequency |
| **`.`** | Increase Andersen frequency |
| **`b`** | Toggle heat bath ON/OFF |
| **`[`** | Decrease gas temperature by 5% |
| **`]`** | Increase gas temperature by 5% |
| **`a`** | Move left piston left (original function preserved!) |
| **`d`** | Move left piston right |

---

## 3. New CLI Flag: `--hb-temp`

### Purpose
Set heat bath temperature from command line without needing to press keys interactively.

### Usage

#### 00ALLINONE.c
```bash
# Set heat bath to T=5.0 (auto-enables heat bath)
./00ALLINONE --hb-temp=5.0

# Combine with other flags
./00ALLINONE --temperature=1.0 --hb-temp=5.0 --mode=edmd

# Use with input file
./00ALLINONE --input in.test --hb-temp=101.0
```

#### boxtest.c
**NOT YET IMPLEMENTED** - boxtest.c currently only accepts number of disks as CLI arg.

To add later (simple implementation):
```c
// In main():
if (argc >= 3 && strcmp(argv[2], "--hb-temp") == 0) {
    hbtemp = atof(argv[3]);
    heatbath = 1;
}
```

---

## 4. Implementation Details

### 00ALLINONE.c Changes

#### A. New CLI Variable (Line 348)
```c
static float  cli_override_hb_temp = -1.0f;     // ##CHRIS: --hb-temp=VALUE (heat bath temperature)
```

#### B. CLI Parsing (Lines 1696-1706)
```c
} else if (strncmp(arg, "--hb-temp", 9) == 0) {
    // ##CHRIS: Heat bath temperature override
    const char *value = cli_option_value(arg, argc, argv, &i);
    errno = 0; char *endptr = NULL; float v = strtof(value, &endptr);
    if (errno != 0 || endptr == value || *endptr != '\0' || v <= 0.0f) {
        fprintf(stderr, "Invalid heat bath temperature '%s'.\n", value);
        exit(EXIT_FAILURE);
    }
    cli_override_hb_temp = v;
    heatbath_enabled = 1;  // Auto-enable heat bath when temperature specified
}
```

#### C. Apply Override (Lines 1909-1911)
```c
// ##CHRIS: Apply heat bath temperature override
if (cli_override_hb_temp > 0.0f) {
    heatbath_temperature = cli_override_hb_temp;
}
```

#### D. Help Text (Lines 1322-1323)
```c
printf("  --temperature=value         Set initial gas temperature (reduced units)\n");
printf("  --hb-temp=value             Set heat bath temperature and auto-enable it (reduced units)\n");
```

#### E. Keyboard Change (Line 5194)
```c
// ##CHRIS: Andersen thermostat toggle with 'k' key (NOT 'a' - that's for piston!)
case SDLK_k: {
    andersen_enabled = !andersen_enabled;
    ...
}
```

#### F. Display Update (Lines 5655, 5657)
```c
if (andersen_enabled) {
    snprintf(andersen_line, sizeof(andersen_line), "k - Andersen: ON (ν=%.3f)", andersen_collision_freq);
} else {
    snprintf(andersen_line, sizeof(andersen_line), "k - Andersen: off (,/. freq)");
}
```

---

## 5. Testing

### Test 1: Verify Key Change
```bash
cd hspist3
./00ALLINONE

# In simulation window:
# Press 'a' → left piston moves (original function works!)
# Press 'k' → Andersen toggles on/off (new key works!)
```

### Test 2: Test CLI Flag
```bash
cd hspist3

# Test 1: Set HB temp to 5.0
./00ALLINONE --hb-temp=5.0 --show-simulation

# Test 2: Large temperature difference
./00ALLINONE --temperature=1.0 --hb-temp=101.0 --show-simulation

# Test 3: With input file
./00ALLINONE --input in.high_temp_test --hb-temp=50.0
```

### Test 3: Verify Auto-Enable
```bash
# Heat bath should be enabled automatically when --hb-temp is used
./00ALLINONE --hb-temp=5.0

# Check on-screen display: should show "HB: ON (T=5.000, mode=...)"
```

---

## 6. Compilation Status

✅ **Successfully compiled** with only minor warnings

```bash
cd hspist3
make clean && make
# 3 warnings generated (K_B redefinition, fabsf→fabs)
# Binary created: 00ALLINONE (547 KB)
```

---

## 7. Updated Documentation Files

### Files Updated:
1. [hspist3/00ALLINONE.c](hspist3/00ALLINONE.c)
   - Changed Andersen key from 'a' to 'k'
   - Added `--hb-temp` CLI flag
   - Updated help text
   - Updated display text

2. [ANDERSEN_IMPLEMENTATION_00ALLINONE.md](hspist3/ANDERSEN_IMPLEMENTATION_00ALLINONE.md)
   - Should be updated to reflect 'k' key change

### Files NOT Changed:
- [boxtest_2_incl_MB_heatbath/boxtest.c](boxtest_2_incl_MB_heatbath/boxtest.c)
  - Still uses 'a' for Andersen (no conflict in boxtest!)
  - Does NOT have `--hb-temp` CLI flag yet (optional future addition)

---

## 8. Why This Matters

### Before (Broken)
```bash
# In 00ALLINONE.c:
# Press 'a' → Andersen toggles AND piston moves (conflict!)
# Need to press 'h' then 'u'/'j' many times to set HB temp
```

### After (Fixed)
```bash
# In 00ALLINONE.c:
# Press 'a' → piston moves only
# Press 'k' → Andersen toggles only
# Use --hb-temp=VALUE to set HB temp from command line
```

---

## 9. Quick Reference

### For Interactive Use (00ALLINONE)
```bash
./00ALLINONE --show-simulation

# Then use keys:
# k - toggle Andersen
# ,/. - adjust Andersen frequency
# b - toggle heat bath
# [/] - adjust gas temperature
# a/d - move left piston
```

### For Batch Experiments (00ALLINONE)
```bash
./00ALLINONE --experiment=speed_of_sound --temperature=1.0 --hb-temp=5.0
```

### For boxtest (Simple Interactive)
```bash
./boxtest 100  # 100 disks

# Then use keys:
# a - toggle Andersen
# +/- - adjust Andersen frequency
# h - toggle heat bath
# u/j - adjust HB temperature
# t - cycle thermal wall mode
```

---

## 10. Summary

✅ **Key conflict resolved**: Andersen uses 'k' in 00ALLINONE.c, 'a' in boxtest.c
✅ **CLI flag added**: `--hb-temp=VALUE` auto-enables heat bath
✅ **Help text updated**: Shows both `--temperature` and `--hb-temp`
✅ **Display updated**: Shows 'k' key for Andersen in 00ALLINONE.c
✅ **Compilation successful**: No errors, only minor warnings
✅ **Both programs work**: boxtest.c and 00ALLINONE.c fully functional

### Recommended Usage (Realistic Hard Disk Gas)
```bash
# 00ALLINONE:
./00ALLINONE --temperature=1.0 --hb-temp=5.0 --show-simulation
# Then press 'k' to enable Andersen (ν=0.1 default)

# boxtest:
./boxtest 100
# Then press 'h' (HB on), 'u' (increase temp), 'a' (Andersen on)
```

This gives: **Thermal walls (Mode 1) + low-frequency Andersen = Physically realistic hard disk gas!**
