// include LiquidCrystal for the LCD
#include <LiquidCrystal.h>
// initialize the class named lcd
LiquidCrystal lcd(6, 7, 8, 9, 10, 11, 12, A1, A2, A3); // LiquidCrystal lcd(rs, enable, d0, d1, d2, d3, d4, d5, d6, d7) // stable


// include DHT for the temperature sensor
#include <dht_nonblocking.h>
DHT_nonblocking dht_sensor(A5, DHT_TYPE_11); // dht_sensor(sensor_pin, sensor_type); the type of the temperature sensor we use is DHT_TYPE_11


// include IRremote for the IR receiver
#include <IRremote.h>      // it is using Timer2, PWM in D3, and LED in D13
IRrecv irrecv(A4);         // irrecv(sensor_pin);  create instance of 'irrecv'; 
decode_results results;    // create instance of 'decode_results'


// pins for the DC motor
#define PIN_ENABLE 5
#define DIRA 2
#define DIRB 4


// motor speed constants
# define HIGH_SPEED 255 // 100% duty
# define MED_SPEED  170  // 66% duty cycle
# define LOW_SPEED  127   // 50% duty cycle
# define NO_SPEED  0     // 0% duty cycle


// global variables
#define DEFAULT_TEMP 25.0       // default desired temperature
byte MENU, pre_MENU;            // used to print respective menu
byte OMD, pre_OMD, prepre_OMD;  // used to save the operation mode
// used to save the operation time
byte OD, OH, OM;
volatile byte OS;
bool ERR1;                      // a flag for error level 1 
bool FLAG_FR;                   // a flag for the first running
unsigned long pre_temp_time;    // used to read temperature in the time interval
float temp_goal, temperature;   // used to save temperature
float humidity;                 // used to save humidity
byte motor_speed;               // used to control the speed of the motor
int timer_goal;                 // used to save timer
bool FLAG_FRS;                  // a flag for the first running setup
bool FLAG_TIMER;                // a flag for the timer
bool skipFlag = false;


void mainmenu() {
  Serial.println("EME-154 Mechatronics");
  Serial.println("Final Project");
  Serial.println("Michael Gunnarson");
  Serial.println("***************************************");
  Serial.println("               Main Menu               ");
  Serial.println("***************************************");
  Serial.println("FUNC/STOP: Idle (Stop) Mode");
  Serial.println("1: Manual Mode");
  Serial.println("2: Automatic Mode");
  Serial.println("3: Setup Mode");
  Serial.println("POWER: Exit Mode");
  Serial.println("");
  delay(50); // wait for printing
}

void manualmenu() {
  Serial.println("***************************************");
  Serial.println("              Manual Menu              ");
  Serial.println("***************************************");
  Serial.println("1: Low Speed");
  Serial.println("2: Medium Speed");
  Serial.println("3: High Speed");
  Serial.println("4: Back to Main Menu");
  Serial.println("");
  delay(50); // wait for printing
}

void setupmenu() {
  Serial.println("**************************************");
  Serial.println("              Setup Menu              ");
  Serial.println("**************************************");
  Serial.println("1: Timer Setup");
  Serial.println("2: Temperature Setup");
  Serial.println("3: Back to Main Menu");
  Serial.println("");
  delay(50); // wait for printing
}

void settimer() {
  // set up Timer
  noInterrupts(); // stop interrupts
  
  //set Timer1 interrupt at 1 Hz = 1 time per second
  TCCR1A = 0; // initialize TCCR1A register
  TCCR1B = 0; // initialize TCCR1B register
  TCNT1  = 0; // initialize counter value to 0
  // set compare match register for 1hz increments
  OCR1A = 15624; // = 16,000,000 / 1 Hz / 1024 prescaler - 1 (must be <65536)
  // turn on CTC mode
  TCCR1B |= (1 << WGM12); // set WGM12 = 1
  // set CS12 and CS10 equal to 1 for 1024 prescaler
  TCCR1B |= (1 << CS12) | (1 << CS10);
  // enable timer compare match A interrupt
  TIMSK1 |= (1 << OCIE1A); // set OCIE1A = 1
  
  interrupts(); // allow interrupts
}

// timer interrupt function with compare match A interrupt in Timer1
ISR(TIMER1_COMPA_vect){
  if (OMD != 0 && !FLAG_TIMER) { // check for idle mode and countdown flag
    OS++;
    if (OS > 59){  // 60 sec in a minute
      OM++;
      OS = 0;
    }
    if (OM > 59){  // 60 minutes in an hour
      OH++;
      OM = 0;
    }
    if (OH > 23){  // 24 hours in a day
      OD++;
      OH = 0;
    }
  }

  
  else if (OMD != 0 && FLAG_TIMER) { // decrement vars
    if (OS ==1 && OM == 0 && OH == 0 && OD == 0) {
      OMD = 0;
      resetTimer();
    }
    else {
      OS--;
    if (OH == 0 && OM == 0 && OS == 0){  // 24 hours in a day
      OD--;
      OH = 24;
    }
    if (OM == 0 && OS == 0){  // 60 minutes in an hour
      OH--;
      OM = 59;
    }
    if (OS == 0){  // 60 sec in a minute
      OM--;
      OS = 59;
    }
  }
  }
  

  // print operation time on the LCD
  lcd.setCursor(0, 1); // lcd.setCursor(col, row)
  lcd.print(OD);
  lcd.print("D");
  lcd.print(OH);
  lcd.print("H");
  lcd.print(OM);
  lcd.print("M");
  lcd.print(OS);
  lcd.print("S  ");

  // put the timer program here to stop the motor when counting down to zero
  
}

// INZ
void INZfunction() {
  // initailize parameters
  MENU = 0;
  pre_MENU = 0;
  OMD = 0;
  pre_OMD = 0;
  prepre_OMD = 0;
  OD = 0;
  OH = 0;
  OM = 0;
  OS = 0;
  ERR1 = false;
  FLAG_FR = true;
  motor_speed = 0;
  temp_goal = DEFAULT_TEMP;
  pre_temp_time = millis();
  FLAG_FRS = true;
  timer_goal = 0;
  FLAG_TIMER = false;
}

// DIG
int DIGfunction() {
  // put the diagnostics program here
  if (ERR1) {
    return 0; // for error
  } else {
    return 1; // successfull diagnostics
  }
}

// ERT
void ERTfunction() {
  // take actions here to treat errors detected
  ERR1 = false;
}

// MSS
void MSSfunction() {
  // scan the IR receiver to get the input
  if (MENU == 0) {
    if (irrecv.decode(&results)) {
      switch (results.value) {
        case 0xFFE21D: // FUNC/STOP
          OMD = 0; // idle mode
        break;
        case 0xFF30CF: // 1
          MENU = 1; // manual mode
        break;
        case 0xFF18E7: // 2
          OMD = 4; // automatic mode
        break;
        case 0xFF7A85: // 3
          MENU = 2; // setup mode
        break;
        case 0xFFA25D: // POWER
          OMD = 7; // turn off system
        break;
        default:
          //Serial.println(results.value); // print the key we press for checking
        break;
      }
      irrecv.resume(); // receive the next value
    }
  } else if (MENU == 1) {
    if (irrecv.decode(&results)) {
      switch (results.value) {
        case 0xFF30CF: // 1
          // Low Speed
          OMD = 1;
        break;
        case 0xFF18E7: // 2
          // Medium Speed
          OMD = 2;
        break;
        case 0xFF7A85: // 3
          // High Speed
          OMD = 3;
        break;
        case 0xFF10EF: // 4
          // back to main menu
          MENU=0;
        break;
        default:
          //Serial.println(results.value); // print the key we press for checking
        break;
      }
      irrecv.resume(); // receive the next value
    }
  } else if (MENU == 2) { // setup mode
    if (irrecv.decode(&results)) {
      switch (results.value) {
        case 0xFF30CF: // 1
          OMD = 5; // setup timer
        break;
        case 0xFF18E7: // 2
          OMD = 6; // setup temp
        break;
        case 0xFF7A85: // 3
          MENU = 0; // main menu
        break;
        default:
          //Serial.println(results.value); // print the key we press for checking
        break;
      }
      irrecv.resume(); // receive the next value
    }
  }

  // get temperature
  if ((millis()-pre_temp_time) > 3000) { // read temperature and humidity every 3 seconds
    if (dht_sensor.measure(&temperature, &humidity) == true) {
      //Serial.print("T = ");
      //Serial.print(temperature);
      //Serial.print(" deg. C, H = ");
      //Serial.print(humidity);
      //Serial.println(" %");
      pre_temp_time = millis();
    }
  }
}

// MCS
void MCSfunction() {
  switch (OMD) {
    case 0:
      IDMfunction(); // Idle (Stop) Mode
    break;
    case 1:
      MOSfunction(); // Manual Mode (Low Speed)
    break;
    case 2:
      MOSfunction(); // Manual Mode (Medium Speed)
    break;
    case 3:
      MOSfunction(); // Manual Mode (High Speed)
    break;
    case 4:
      ACMfunction(); // Automatic Mode
    break;
    case 5:
      MSDfunction(); // Setup Mode (Timer Setup)
    break;
    case 6:
      MSDfunction(); // Setup Mode (Temperature Setup)
    break;
    case 7:
      // System Turned Off;
      IDMfunction(); // Idle (Stop) Mode
    break;
    default:
    break;
  }
}

// Idle (Stop) Mode
void IDMfunction() {
  // the idle (stop) mode program
  // put stop motor and reset timer program here
  NSMfunction();
}

// MOS
void MOSfunction() {
  switch (OMD) {
    case 1:
      // Low Speed
      LSMfunction();
    break;
    case 2:
      // Medium Speed
      MSMfunction();
    break;
    case 3:
      // High Speed
      HSMfunction();
    break;
    default:
    break;
  }
}

void setMotorCW(){
  digitalWrite(DIRA,LOW); //one way
  digitalWrite(DIRB,HIGH);
}

// Low Speed Mode
void LSMfunction() {
  // the low speed mode program
  // put low speed program here
  setMotorCW();
  analogWrite(PIN_ENABLE, LOW_SPEED);
}

// Medium Speed Mode
void MSMfunction() {
  // the medium speed mode program
  // put medium speed program here
  setMotorCW();
  analogWrite(PIN_ENABLE, MED_SPEED);
}

// High Speed Mode
void HSMfunction() {
  // the high speed mode program
  // put high speed program here
  setMotorCW();
  analogWrite(PIN_ENABLE, HIGH_SPEED);
}

// No Speed Mode
void NSMfunction() {
  // the no speed mode program
  // put no speed program here
  setMotorCW();
  analogWrite(PIN_ENABLE, NO_SPEED);
}

// Automatic Mode
void ACMfunction() {
  // the automatic mode program
  // put automatic speed program here
  float error = temperature - temp_goal;
  //Serial.println(error);
  if(error >= 1.0) {
    HSMfunction();
  }
  else if (error >= 0.5) {
    MSMfunction();
  }
  else if (error >= 0.0) {
    LSMfunction();
  }
  else {
    NSMfunction();
  }
}

// MSD
void MSDfunction() {
  switch (OMD) {
    case 5:
      TISMfunction(); // Timer Setup
    break;
    case 6:
      TPSMfunction(); // Temperature Setup
    break;
    default:
    break;
  }
}

// Timer Setup Mode
void TISMfunction() {
  // the timer setup mode program
  if (!FLAG_FRS) {
    if (prepre_OMD != 0) {
      // put the timer setup process here
      Serial.println("Insert the time you want the motor to run in seconds into the serial terminal.");
      Serial.println("Terminate input with enter button");
      timer_goal = getUserInput();
      setTimer(timer_goal);
      Serial.print(timer_goal); Serial.println(" seconds");
      prepre_OMD = 4;
      OMD = 4; // start automatic mode
      
    } else {
      //Serial.println("The system is in the idle (stop) mode already!");
      //Serial.println("Timer Setup Failed!");
      //Serial.println("");
      Serial.println("Insert the time you want the motor to run in seconds into the serial terminal.");
      Serial.println("Terminate input with enter button");
      timer_goal = getUserInput();
      setTimer(timer_goal);
      Serial.print(timer_goal); Serial.println(" seconds");
      prepre_OMD = 4;
      OMD = 4; // start automatic mode
    }
    FLAG_FRS = !FLAG_FRS;
    OMD = prepre_OMD;
  } else {
    FLAG_FRS = !FLAG_FRS;
  }
}

int getUserInput(){
  int number;
  while(1) {
    if(Serial.available() > 0) {
      String str = Serial.readStringUntil('\n');
      number = str.toInt();
      //Serial.print(number);
      break;
    }
  }
return number;
}

void resetTimer(){
  FLAG_TIMER = false;
  OD = 0;
  OH = 0;
  OM = 0;
  OS = 0;
}

void setTimer(int seconds){
  FLAG_TIMER = true;

  long int secRem = seconds;
  int minRem = 0;
  int hourRem = 0;
  int minCount = 0;
  int hourCount = 0;
  int dayCount = 0;
  while(secRem-60>=0){
    minCount += 1;
    secRem = secRem-60;
  }
  minRem = minCount;
  while(minRem-60>=0) {
    hourCount += 1;
    minRem = minRem - 60;
  }
  hourRem = hourCount;
  while(hourRem-24>=0){
    dayCount += 1;
    hourRem = hourRem - 24;
  }
  
  int day = dayCount;
  int hour = hourRem;
  int minute = minRem;
  int second = secRem;
   
  OD = day;
  OH = hour;
  OM = minute;
  OS = second;
}

// Temperature Setup Mode
void TPSMfunction() {
  // the temperature setup mode program
  if (!FLAG_FRS) {
    // put the temperature setup process here
      Serial.println("Insert the temperature you want in Celsius");
      Serial.println("Terminate input with enter button");
      temp_goal = getUserInput();
    FLAG_FRS = !FLAG_FRS;
    OMD = prepre_OMD;
  } else {
    FLAG_FRS = !FLAG_FRS;
  }
}



char *getAcceptedMode() {
  switch (OMD) {
    case 0:
      return "Idle (Stop) Mode Accepted";
    case 1:
      return "Low Speed Mode Accepted";
    case 2:
      return "Medium Speed Mode Accepted";
    case 3:
      return "High Speed Mode Accepted";
    case 4:
      return "Automatic Mode Accepted";
    case 5:
      return "Timer Setup Mode Accepted";
    case 6:
      return "Temperature Setup Mode Accepted";
    case 7:
      return "Exit Command Accepted";
    default:
      return "Stop";
  }
}

char *getCurrentMode() {
  switch (OMD) {
    case 0:
      return "Idle (Stop)";
    case 1:
      return "Low";
    case 2:
      return "Medium";
    case 3:
      return "High";
    case 4:
      return "Automatic";
    case 5:
      return "Timer Setup";
    case 6:
      return "Temperature Setup";
    case 7:
      return "System Off";
    default:
      return "Stop";
  }
}

// OCS
void OCSfunction() {
  // print on the serial monitor
  // print menu
  if (pre_MENU != MENU || FLAG_FR == true) {
    if (MENU == 0) {
      mainmenu();
    } else if (MENU == 1) {
      manualmenu();
    } else if (MENU == 2) {
      setupmenu();
    }
    pre_MENU = MENU;
  }

  if (pre_OMD != OMD || FLAG_FR == true) {
    Serial.println(getAcceptedMode());
    Serial.print("Current Mode: ");
    Serial.println(getCurrentMode());
    Serial.println("");
    prepre_OMD = pre_OMD;
    pre_OMD = OMD;
    FLAG_FR = false;
    delay(50); // wait for printing
  }

  // print on the LCD
  // print temperature and humidity
  lcd.setCursor(0, 0); // lcd.setCursor(col, row)
  lcd.print("T:");
  lcd.print(int(temperature));
  lcd.print(" H:");
  lcd.print(int(humidity));
  lcd.print(" DT:");
  lcd.print(int(temp_goal));
  
  // print timer
  lcd.setCursor(10, 1); // lcd.setCursor(col, row)
  lcd.print("TR:");
  lcd.print(timer_goal);
  lcd.print("    ");
  
  // print modes
  switch (OMD) {
    case 0:
      lcd.setCursor(15, 0); // lcd.setCursor(col, row)
      lcd.print("I");
    break;
    case 1:
      lcd.setCursor(15, 0); // lcd.setCursor(col, row)
      lcd.print("M");
    break;
    case 2:
      lcd.setCursor(15, 0); // lcd.setCursor(col, row)
      lcd.print("M");
    break;
    case 3:
      lcd.setCursor(15, 0); // lcd.setCursor(col, row)
      lcd.print("M");
    break;
    case 4:
      lcd.setCursor(15, 0); // lcd.setCursor(col, row)
      lcd.print("A");
    break;
    case 5:
      lcd.setCursor(15, 0); // lcd.setCursor(col, row)
      lcd.print("S");
    break;
    case 6:
      lcd.setCursor(15, 0); // lcd.setCursor(col, row)
      lcd.print("S");
    break;
    case 7:
      noInterrupts();
      lcd.setCursor(0, 0); // lcd.setCursor(col, row)
      lcd.print(getCurrentMode());
      lcd.print("!     ");
      lcd.setCursor(0, 1); // lcd.setCursor(col, row)
      lcd.print("Exit!           ");
    break;
    default:
    break;
  }
}

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  lcd.begin(16, 2); // lcd.begin(cols, rows)
  lcd.clear(); // clear the LCD screen
  irrecv.enableIRIn(); // enable the IR receiver
  // set functions of pins for the motor
  pinMode(PIN_ENABLE, OUTPUT);
  pinMode(DIRA, OUTPUT);
  pinMode(DIRB, OUTPUT);
  digitalWrite(DIRA, LOW);
  digitalWrite(DIRB, LOW);
  // set the timer
  settimer();
  // the initialization
  INZfunction();
}

void loop() {
  // put your main code here, to run repeatedly:
  // the main FTS
  if (DIGfunction()) {
    // if diagnostics succeed
    MSSfunction();
    MCSfunction();
    OCSfunction();
  } else {
    // go to ERT if diagnostics fail
    ERTfunction();
    OCSfunction();
  }
  // if press POWER, stop the system
  if (OMD == 7) {
    while (1) {};
  }
}
