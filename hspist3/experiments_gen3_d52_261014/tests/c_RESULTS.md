# The return-map panel's tests (c), the binary in tests/c_H_e030_closed.cmd

## H, eta 0.30 (N 400): gen3's event hash at exactly steps 600, 1200 and 1800

| run | step 600: hash, t, events | step 1200 | step 1800 | panel notes | run record clean |
|---|---|---|---|---|---|
| closed | 0962febcf90cda92, 240.00000357627869, 10163 | 8d360cfc76abd82b, 480.00000715255737, 20538 | f92dfdfd34f0504f, 720.00001072883606, 30918 | 0:  | 1 |
| open | 0962febcf90cda92, 240.00000357627869, 10163 | 8d360cfc76abd82b, 480.00000715255737, 20538 | f92dfdfd34f0504f, 720.00001072883606, 30918 | 1: return map: angle line of centres, walls shown, side both, K 20000 | 1 |
| toggled | 0962febcf90cda92, 240.00000357627869, 10163 | 8d360cfc76abd82b, 480.00000715255737, 20538 | f92dfdfd34f0504f, 720.00001072883606, 30918 | 7: return map: angle line of centres, walls shown, side both, K 20000; I: return map closed; return map: angle line of centres, walls shown, side both, K | 1 |
| keys | 0962febcf90cda92, 240.00000357627869, 10163 | 8d360cfc76abd82b, 480.00000715255737, 20538 | f92dfdfd34f0504f, 720.00001072883606, 30918 | 12: return map: angle line of centres, walls shown, side both, K 20000; return map: angle velocity after, walls shown, side both, K 20000; return map: ang | 1 |

verdict: EQUAL (the four runs at the three steps)

## H, eta 0.69 (N 400): gen3's event hash at exactly steps 600, 1200 and 1800

| run | step 600: hash, t, events | step 1200 | step 1800 | panel notes | run record clean |
|---|---|---|---|---|---|
| closed | 0047e2d2bca75315, 240.00000357627869, 62490 | dcedee2ec1966e1e, 480.00000715255737, 123507 | 05d74569dce5416e, 720.00001072883606, 185391 | 0:  | 1 |
| open | 0047e2d2bca75315, 240.00000357627869, 62490 | dcedee2ec1966e1e, 480.00000715255737, 123507 | 05d74569dce5416e, 720.00001072883606, 185391 | 1: return map: angle line of centres, walls shown, side both, K 20000 | 1 |
| toggled | 0047e2d2bca75315, 240.00000357627869, 62490 | dcedee2ec1966e1e, 480.00000715255737, 123507 | 05d74569dce5416e, 720.00001072883606, 185391 | 7: return map: angle line of centres, walls shown, side both, K 20000; I: return map closed; return map: angle line of centres, walls shown, side both, K | 1 |
| keys | 0047e2d2bca75315, 240.00000357627869, 62490 | dcedee2ec1966e1e, 480.00000715255737, 123507 | 05d74569dce5416e, 720.00001072883606, 185391 | 12: return map: angle line of centres, walls shown, side both, K 20000; return map: angle velocity after, walls shown, side both, K 20000; return map: ang | 1 |

verdict: EQUAL (the four runs at the three steps)

## F, the frames (PNG in tests/)

- c_F_e030_centres: 2000 x 1200 px window -> tests/c_F_e030_centres.png and tests/c_F_e030_centres_panel.png (the panel)
- c_F_e069_centres: 2000 x 1200 px window -> tests/c_F_e069_centres.png and tests/c_F_e069_centres_panel.png (the panel)
- c_F_e069_velocity: 2000 x 1200 px window -> tests/c_F_e069_velocity.png and tests/c_F_e069_velocity_panel.png (the panel)
