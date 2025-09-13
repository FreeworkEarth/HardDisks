00_EDMD_TIME/
├─ core/
│  ├─ state.h                // shared Particle/Params/State structs
│  ├─ edmd_core.c/.h         // event-driven engine (CPU)
│  ├─ step_core.c/.h         // your old time-step logic, split out
│  ├─ protocol.c/.h          // scripted piston/wall protocols
│  └─ io_log.c/.h            // CSV/TSV logging helpers
├─ ui/
│  ├─ viz.c/.h               // SDL2 renderer (toggle on/off)
│  └─ controls.c/.h          // keyboard → control signals
├─ app/
│  ├─ main.c                 // CLI, engine switch, run loop
│  └─ config_cli.c/.h        // parse flags → Params/Protocol config
├─ third_party/              // (optional) imgui, stb, etc.
└─ build scripts (Makefile or CMake)