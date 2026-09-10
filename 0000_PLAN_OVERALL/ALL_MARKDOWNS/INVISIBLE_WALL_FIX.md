# Invisible Wall Fix: W Key Now Syncs with EDMD

## The Problem

User reported: *"the wall somehow although i switch it off pressing W and it dissappears, the wall is still there for the particles and they intereact with the wall in edmd mode"*

**Root cause**: The 'W' key toggled `wall_enabled` for **rendering only**, but did NOT synchronize with EDMD's internal `has_divider` parameter!

### What Was Happening

1. User presses 'W' to disable wall
2. `wall_enabled = 0` → Wall becomes invisible in graphics
3. **EDMD still has `has_divider = 1`** → Particles still collide with invisible wall!
4. EDMD's collision schedule was never updated

**Result**: Ghost wall - visually gone, but physically still blocking particles.

---

## The Fix

### File Modified: [hspist3/00ALLINONE.c](hspist3/00ALLINONE.c)

**Lines 5136-5157** - Updated 'W' key handler with EDMD synchronization

### Before (Broken)

```c
case SDLK_w:  // Toggle wall with 'w' or 'W'
    wall_enabled = !wall_enabled;
    printf("Wall %s\n", wall_enabled ? "Enabled" : "Disabled");
    break;
```

**Problem**: Only toggles rendering flag, EDMD unaware!

### After (Fixed)

```c
case SDLK_w:  // ##CHRIS: Toggle wall with 'w' or 'W' (syncs with EDMD)
    wall_enabled = !wall_enabled;

    // ##CHRIS: Sync wall state with EDMD's divider parameter
    if (g_edmd) {
        // Get current divider parameters
        const EDMD_Params* ep = edmd_params(g_edmd);

        // Update divider state (enabled/disabled)
        edmd_config_divider(g_edmd,
                           wall_enabled ? 1 : 0,  // enabled flag
                           ep->divider_x,          // keep current position
                           ep->divider_thickness); // keep current thickness

        // Rebuild collision schedule so EDMD knows about the change
        edmd_reschedule_all(g_edmd);

        printf("Wall %s (EDMD divider synced)\n", wall_enabled ? "Enabled" : "Disabled");
    } else {
        printf("Wall %s\n", wall_enabled ? "Enabled" : "Disabled");
    }
    break;
```

**Solution**:
1. Update EDMD's `has_divider` parameter using `edmd_config_divider()`
2. Rebuild collision schedule with `edmd_reschedule_all()`
3. Print confirmation message showing EDMD was synced

---

## How It Works

### Step 1: Get Current EDMD Parameters
```c
const EDMD_Params* ep = edmd_params(g_edmd);
```
Read current divider position and thickness (we want to keep these unchanged)

### Step 2: Update Divider State
```c
edmd_config_divider(g_edmd,
                   wall_enabled ? 1 : 0,  // NEW: enabled/disabled flag
                   ep->divider_x,          // KEEP: current position
                   ep->divider_thickness); // KEEP: current thickness
```
Only change the `has_divider` flag, preserve position/thickness

### Step 3: Rebuild Collision Schedule
```c
edmd_reschedule_all(g_edmd);
```

**Why needed?** EDMD maintains a priority queue of upcoming collision events. When the divider is removed:
- All particle-divider collision events must be removed from queue
- Particles that were about to hit the divider can now pass through
- New particle-particle collisions may be scheduled

**What `edmd_reschedule_all()` does**:
1. Clears the event priority queue
2. Rebuilds spatial grid for collision detection
3. Reschedules all particle-particle collisions
4. Reschedules all wall collisions (including divider, if enabled)

From [edmd.c:306-319](hspist3/edmd_core/edmd.c):
```c
static void reschedule_all_internal(EDMD* S){
    heap_clear(&S->pq);           // Clear event queue
    rebuild_grid(S);              // Rebuild spatial grid
    for(int i=0;i<S->prm.N;i++){
        schedule_walls(S,i);      // Schedule wall collisions (includes divider)
        schedule_divider(S,i);    // Schedule divider collisions if has_divider=1
        schedule_pistons(S,i);    // Schedule piston collisions if enabled
        schedule_cells(S,i);      // Schedule particle-particle collisions
    }
}
```

---

## Testing

### Test 1: EDMD Mode - Toggle Wall Off
```bash
cd hspist3
./00ALLINONE --mode=edmd --particles=100 --l0=20 --height=10 \
  --wall-mass-factor=200 --temperature=137000 --kbt1 --hb-temp=500000 \
  --show-simulation

# During simulation:
# 1. Observe particles bouncing off wall
# 2. Press 'W' to disable wall
# 3. Wall disappears AND particles pass through! ✅
# 4. Press 'W' again to re-enable
# 5. Wall reappears AND particles bounce again! ✅
```

**Expected output**:
```
Wall Disabled (EDMD divider synced)  ← Confirmation message
Wall Enabled (EDMD divider synced)   ← Confirmation message
```

### Test 2: Non-EDMD Mode - Still Works
```bash
./00ALLINONE --mode=ccd --particles=100 --show-simulation

# Press 'W' to toggle wall
# Should still work (no EDMD to sync)
```

**Expected output**:
```
Wall Disabled  ← Simple message (no EDMD)
Wall Enabled
```

### Test 3: Your Original Command
```bash
./00ALLINONE --mode=edmd --particles=100 --l0=20 --height=10 \
  --wall-mass-factor=200 --temperature=137000 --kbt1 --hb-temp=500000

# Press 'W' → Wall should truly disappear for particles now!
```

---

## Why the Pattern Was Copied from Heat Bath Toggle

The fix follows the **same pattern** as the heat bath toggle ('b' key) at lines 5189-5208:

| Heat Bath Toggle | Wall Toggle |
|------------------|-------------|
| `heatbath_enabled = !heatbath_enabled;` | `wall_enabled = !wall_enabled;` |
| `if (g_edmd) {` | `if (g_edmd) {` |
| `((EDMD_Params*)edmd_params(g_edmd))->heatbath_enabled = ...` | `edmd_config_divider(g_edmd, wall_enabled ? 1 : 0, ...)` |
| No reschedule needed (changes velocity on next collision) | `edmd_reschedule_all(g_edmd);` (must remove/add divider events) |

**Key difference**: Heat bath doesn't need rescheduling because it only changes **how** wall collisions are resolved, not **which** collisions exist. Divider toggle changes **which** collisions exist, so reschedule is required.

---

## Implementation Notes

### Function Signatures Used

From [hspist3/edmd_core/edmd.h](hspist3/edmd_core/edmd.h):

```c
/* Line 63: Get EDMD parameters (read-only pointer) */
const EDMD_Params* edmd_params(const EDMD* S);

/* Line 77: Rebuild collision schedule after bulk changes */
void edmd_reschedule_all(EDMD* S);

/* Line 84: Configure divider slab */
void edmd_config_divider(EDMD* S, int enabled, double cx, double thickness);
```

### EDMD_Params Fields (from edmd.h:22-24)

```c
int    has_divider;       /* 1 to enable, 0 to disable */
double divider_x;         /* center x position (0..boxW) */
double divider_thickness; /* slab thickness (>=0) */
```

---

## Compilation Status

✅ **Successfully compiled** (no new errors or warnings)

```bash
cd hspist3
make clean && make
# 3 warnings generated (pre-existing: K_B redefinition, fabsf→fabs)
# Binary created: 00ALLINONE (547 KB)
```

---

## Related Documentation

See also:
- [EDMD_HEATBATH_FIX.md](EDMD_HEATBATH_FIX.md) - Heat bath mode synchronization with EDMD
- [KEY_FIX_AND_CLI_UPDATE.md](KEY_FIX_AND_CLI_UPDATE.md) - Andersen key fix and --hb-temp CLI flag

---

## Summary

✅ **Wall toggle now fully functional in EDMD mode**
✅ **'W' key syncs both rendering AND physics**
✅ **EDMD collision schedule automatically rebuilt**
✅ **Particles can now pass through disabled walls**
✅ **Non-EDMD modes still work (backwards compatible)**

**The ghost wall is gone!** 👻🚫

### Before Fix
```
Press 'W' → Wall invisible, particles still bounce (ghost wall!)
```

### After Fix
```
Press 'W' → Wall invisible AND passable (truly disabled!)
```

**Try your command again - wall toggling should work properly now!** 🎉

```bash
./00ALLINONE --mode=edmd --particles=100 --l0=20 --height=10 \
  --wall-mass-factor=200 --temperature=137000 --kbt1 --hb-temp=500000
```

Press 'W' and watch particles pass through the space where the wall was!
