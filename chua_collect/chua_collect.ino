// =============================================================
//  Chua Circuit Data Logger — R4 Minima, Memory-Safe Buffered
//  Channels: A0 = Vc1 (x), A1 = Vc2 (y), A2 = Il sense (z)
//
//  Memory layout (R4 Minima: 32 KB SRAM total):
//   - Buffer: 3 x uint16_t per sample = 6 bytes/sample
//     1000 samples = 6 000 bytes (~6 KB)
//   - Serial TX buffer (BSP): ~4 KB
//   - Stack + BSP globals:    ~8 KB
//   - Remaining headroom:     ~14 KB  <-- safe
//
//  Timestamps are reconstructed from (sample_index * SAMPLE_INTERVAL_US)
//  rather than stored, saving 4 bytes per sample.
//
//  Serial commands (Serial Monitor, line ending = Newline):
//    1-4  select regime (Fixed Point / Limit Cycle / Period-Doubled / Double Scroll)
//    s    start capture into RAM buffer, dumps CSV when done
//    q    abort capture in progress
//    ?    print menu
// =============================================================

#include <Arduino.h>

// ---- Pin assignments ----------------------------------------
const int PIN_X = A0;   // V_C1
const int PIN_Y = A1;   // V_C2
const int PIN_Z = A2;   // I_L sense (gyrator op-amp output)

// ---- Voltage scaling ----------------------------------------
const float VREF    = 5.0f;
const float ADC_MAX = 4095.0f;   // 12-bit

// Bipolar offset: set to 2.5 if you DC-biased Chua signals into 0-5 V.
// Set to 0.0 if signals are already centred at 0 V (rail-to-rail op-amp).
const float OFFSET_V = 2.5f;

// ---- Sampling interval --------------------------------------
// Target ~20 samples per oscillation cycle.
// Formula: SAMPLE_INTERVAL_US = 1 000 000 / (circuit_freq_Hz * 20)
//
//   Circuit ~100 Hz  ->  500  us  (Fs = 2   kHz)
//   Circuit ~500 Hz  ->  100  us  (Fs = 10  kHz)
//   Circuit ~2  kHz  ->   25  us  (Fs = 40  kHz)  <- near R4 floor
//
// R4 floor: 3 x analogRead at 12-bit = ~15-20 us minimum.
// Start at 500 and halve until the phase portrait looks like clean orbits.
const unsigned long SAMPLE_INTERVAL_US = 30;

// ---- Buffer size --------------------------------------------
// 6 bytes/sample (3 x uint16_t).  Keep total <= ~6 KB to leave
// plenty of room for Serial buffers and stack.
// 1000 samples x 6 B = 6 000 B.  Safe on 32 KB R4 Minima.
// If you need more samples, reduce to 800 and re-test.
const int MAX_SAMPLES = 1000;

// Raw ADC values only — timestamps reconstructed on dump.
// Using a plain parallel arrays (not a struct) avoids any alignment padding.
static uint16_t bufX[MAX_SAMPLES];
static uint16_t bufY[MAX_SAMPLES];
static uint16_t bufZ[MAX_SAMPLES];

// Capture start time (micros) — one value, not per-sample
static unsigned long captureStartUs = 0;

// ---- State --------------------------------------------------
volatile bool running  = false;
int           regime   = 0;
int           nSamples = 0;

const char* REGIME_NAMES[] = {
    "",
    "fixed_point",
    "limit_cycle",
    "period_doubled",
    "double_scroll"
};

// ---- ADC setup ----------------------------------------------
void setupFastADC() {
    // Native 12-bit on RA4M1: ~5 us per analogRead, no extra config needed.
    analogReadResolution(12);
}

// ---- Menu ---------------------------------------------------
void printMenu() {
    Serial.println(F("# ================================================"));
    Serial.println(F("# Chua Logger — R4 Minima (buffered, memory-safe)"));
    Serial.println(F("#   1  Regime 1: Fixed Point"));
    Serial.println(F("#   2  Regime 2: Limit Cycle"));
    Serial.println(F("#   3  Regime 3: Period-Doubled"));
    Serial.println(F("#   4  Regime 4: Double Scroll (chaos)"));
    Serial.println(F("#   s  Start capture (silent during sample, dumps after)"));
    Serial.println(F("#   q  Abort"));
    Serial.println(F("#   ?  This menu"));
    Serial.print  (F("# SAMPLE_INTERVAL_US = ")); Serial.println(SAMPLE_INTERVAL_US);
    Serial.print  (F("# MAX_SAMPLES        = ")); Serial.println(MAX_SAMPLES);
    Serial.print  (F("# Capture duration   = "));
    Serial.print  (SAMPLE_INTERVAL_US * (unsigned long)MAX_SAMPLES / 1000UL);
    Serial.println(F(" ms"));
    Serial.println(F("# Buffer RAM used    = "));
    Serial.print  (MAX_SAMPLES * 6);
    Serial.println(F(" bytes"));
    Serial.println(F("# ================================================"));
    if (regime > 0) {
        Serial.print(F("# Selected: "));
        Serial.print(regime);
        Serial.print(F(" — "));
        Serial.println(REGIME_NAMES[regime]);
        Serial.println(F("# Type s to start."));
    } else {
        Serial.println(F("# No regime selected — type 1, 2, 3, or 4."));
    }
}

// ---- Capture (blocking, no Serial I/O during sampling) ------
void runCapture() {
    Serial.println(F("# Capturing — Serial paused during sampling..."));
    Serial.flush();

    nSamples       = 0;
    captureStartUs = micros();
    unsigned long t_next = captureStartUs;

    while (nSamples < MAX_SAMPLES) {
        // Busy-wait for precise interval (more accurate than delayMicroseconds)
        while ((long)(micros() - t_next) < 0) { /* spin */ }
        t_next += SAMPLE_INTERVAL_US;

        bufX[nSamples] = (uint16_t)analogRead(PIN_X);
        bufY[nSamples] = (uint16_t)analogRead(PIN_Y);
        bufZ[nSamples] = (uint16_t)analogRead(PIN_Z);
        nSamples++;

        // Poll for abort every 128 samples — Serial.available() is very fast
        if ((nSamples & 0x7F) == 0 && Serial.available()) {
            char c = Serial.peek();
            if (c == 'q' || c == 'Q') {
                Serial.read();
                running = false;
                Serial.println(F("# Aborted."));
                return;
            }
        }
    }

    running = false;
}

// ---- Convert raw ADC -> voltage and dump as CSV -------------
void dumpBuffer() {
    const float scale = VREF / ADC_MAX;

    Serial.println(F("# --- DATA START ---"));
    Serial.print  (F("# regime: "));            Serial.println(REGIME_NAMES[regime]);
    Serial.print  (F("# regime_id: "));         Serial.println(regime);
    Serial.print  (F("# n_samples: "));         Serial.println(nSamples);
    Serial.print  (F("# sample_interval_us: ")); Serial.println(SAMPLE_INTERVAL_US);
    Serial.print  (F("# fs_hz: "));
    Serial.println(1000000.0f / (float)SAMPLE_INTERVAL_US, 2);
    Serial.print  (F("# capture_start_us: "));  Serial.println(captureStartUs);
    Serial.println(F("time_ms,x_V,y_V,z_V"));

    char row[52];
    for (int i = 0; i < nSamples; i++) {
        // Reconstruct timestamp from start + index (no per-sample storage needed)
        float t_ms = (captureStartUs + (unsigned long)i * SAMPLE_INTERVAL_US) / 1000.0f;
        float vx   = bufX[i] * scale - OFFSET_V;
        float vy   = bufY[i] * scale - OFFSET_V;
        float vz   = bufZ[i] * scale - OFFSET_V;

        snprintf(row, sizeof(row), "%.3f,%.4f,%.4f,%.4f", t_ms, vx, vy, vz);
        Serial.println(row);
    }

    Serial.println(F("# --- DATA END ---"));
    Serial.print  (F("# Samples dumped: ")); Serial.println(nSamples);
    Serial.println(F("# Next: adjust circuit, type 1-4 then s."));
    Serial.println();
}

// ---- Setup --------------------------------------------------
void setup() {
    // R4 USB-CDC: baud value is ignored by the USB stack (runs at USB 2.0
    // full speed ~1 MB/s effective). Set high so the host driver doesn't
    // insert artificial throttling.
    Serial.begin(2000000);
    while (!Serial);

    setupFastADC();
    running = false;
    printMenu();
}

// ---- Loop ---------------------------------------------------
void loop() {
    if (!Serial.available()) return;

    String line = Serial.readStringUntil('\n');
    line.trim();

    if (line.length() == 0) {
        if (!running) printMenu();
        return;
    }

    char cmd = line.charAt(0);

    if (cmd >= '1' && cmd <= '4') {
        if (running) {
            Serial.println(F("# Abort first (q)."));
        } else {
            regime = cmd - '0';
            Serial.print(F("# Regime -> "));
            Serial.print(regime);
            Serial.print(F(": "));
            Serial.println(REGIME_NAMES[regime]);
            Serial.println(F("# Adjust circuit, then type s."));
        }
    }
    else if (cmd == 's' || cmd == 'S') {
        if (regime == 0) {
            Serial.println(F("# Select a regime first (type 1, 2, 3, or 4)."));
        } else if (!running) {
            running = true;
            runCapture();
            if (nSamples > 0) dumpBuffer();
            printMenu();
        } else {
            Serial.println(F("# Already running."));
        }
    }
    else if (cmd == 'q' || cmd == 'Q') {
        running = false;
        Serial.println(F("# Stopped."));
        printMenu();
    }
    else if (cmd == '?') {
        printMenu();
    }
    else {
        Serial.print(F("# Unknown: '"));
        Serial.print(cmd);
        Serial.println(F("' — type ? for help."));
    }
}
