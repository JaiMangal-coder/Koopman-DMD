// =============================================================
//  Chua Circuit Data Logger
//  Channels: A0 = Vc1 (x), A1 = Vc2 (y), A2 = Il sense (z)
//  Output:   CSV over Serial — timestamp_ms, x_V, y_V, z_V
// =============================================================

// ---- Configuration ------------------------------------------
const int PIN_X = A0;           // V_C1
const int PIN_Y = A1;           // V_C2
const int PIN_Z = A2;           // I_L measurement (output of second op amp in gyrator)

const float VREF      = 5.0;    // Arduino reference voltage (5.0 or 3.3)
const int   ADC_BITS  = 1023;   // 10-bit ADC: 0-1023

// Sampling
const unsigned long SAMPLE_INTERVAL_US = 1000; // 1 ms = 1 kHz
                                                // increase for slow circuits
                                                // e.g. 10000 = 100 Hz

// Offset correction (set these if you used a +2.5V bias to shift
// bipolar Chua voltages into 0-5V range for the Arduino ADC)
const float OFFSET_V  = 2.5;    // subtract this to recover signed voltage
                                 // set to 0.0 if your signal is already 0-5V

// How many samples to collect before stopping (0 = run forever)
const unsigned long MAX_SAMPLES = 5000;

// ---- Globals ------------------------------------------------
unsigned long sampleCount = 0;
unsigned long lastSampleTime = 0;
bool running = true;

// ---- Setup --------------------------------------------------
void setup() {
  Serial.begin(115200);
  while (!Serial);              // wait for Serial Monitor to open

  analogReference(DEFAULT);    // use DEFAULT (5V) or INTERNAL (1.1V)
                                // INTERNAL gives better resolution for
                                // small signals but clips above 1.1V

  // Print header
  Serial.println("# Chua Circuit Data Log");
  Serial.println("# Columns: time_ms, x_V (Vc1), y_V (Vc2), z_V (Il)");
  Serial.println("# Send 's' to start, 'q' to stop");
  Serial.println("time_ms,x_V,y_V,z_V");

  running = false;              // wait for start command
}

// ---- Helpers ------------------------------------------------

// Read ADC and convert to voltage, subtract offset for bipolar signals
float readVoltage(int pin) {
  int raw = analogRead(pin);
  float v = (raw / (float)ADC_BITS) * VREF;
  return v - OFFSET_V;
}

// Oversample for better resolution (average N readings)
float readVoltageOS(int pin, int n) {
  long sum = 0;
  for (int i = 0; i < n; i++) sum += analogRead(pin);
  float v = ((sum / (float)n) / ADC_BITS) * VREF;
  return v - OFFSET_V;
}

// ---- Loop ---------------------------------------------------
void loop() {

  // Check for serial commands
  if (Serial.available()) {
    char cmd = Serial.read();
    if (cmd == 's' || cmd == 'S') {
      running = true;
      sampleCount = 0;
      lastSampleTime = micros();
      Serial.println("# --- START ---");
    }
    if (cmd == 'q' || cmd == 'Q') {
      running = false;
      Serial.println("# --- STOP ---");
    }
  }

  if (!running) return;

  // Non-blocking sample timing
  unsigned long now = micros();
  if (now - lastSampleTime < SAMPLE_INTERVAL_US) return;
  lastSampleTime += SAMPLE_INTERVAL_US;

  // Read three channels
  // Use readVoltageOS(pin, 4) instead for 4x oversampling if noise is bad
  float x = readVoltage(PIN_X);
  float y = readVoltage(PIN_Y);
  float z = readVoltage(PIN_Z);

  // Timestamp in milliseconds
  float t_ms = now / 1000.0;

  // Print CSV row
  Serial.print(t_ms, 3);
  Serial.print(',');
  Serial.print(x, 4);
  Serial.print(',');
  Serial.print(y, 4);
  Serial.print(',');
  Serial.println(z, 4);

  sampleCount++;

  // Stop after MAX_SAMPLES if limit set
  if (MAX_SAMPLES > 0 && sampleCount >= MAX_SAMPLES) {
    running = false;
    Serial.println("# --- DONE ---");
    Serial.print("# Collected ");
    Serial.print(sampleCount);
    Serial.println(" samples");
  }
}
