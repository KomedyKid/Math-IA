import math
import tkinter as tk
from tkinter import messagebox
from datetime import datetime, timedelta
import requests
from timezonefinder import TimezoneFinder
import pytz

def calculate_julian_day(Y, M, D, hour=0, minute=0, second=0):
    if M <= 2:
        Y -= 1
        M += 12

    A = math.floor(Y / 100)
    B = 2 - A + math.floor(A / 4)

    # Calculate the Julian Day Number (JD)
    JD = (math.floor(365.25 * (Y + 4716)) +
          math.floor(30.6001 * (M + 1)) +
          D + B - 1524.5)

    # Convert time to a fraction of a day and add to JD
    day_fraction = (hour + minute / 60 + second / 3600) / 24.0
    JD += day_fraction

    return JD

def solar_declination(JD):
    # JD for January 1, 2000, at 12:00 PM (noon)
    JD_2000 = 2451545.0
    days_since_2000 = JD - JD_2000

    # Number of days in the year
    N = days_since_2000 % 365.24

    # Calculate solar declination in radians
    declination_radians = math.asin(math.sin(math.radians(-23.44)) *
                                    math.cos(math.radians((360 / 365.24) * (N + 10))))

    # Convert to degrees
    declination_degrees = math.degrees(declination_radians)

    return declination_degrees

def calculate_right_ascension(JD):
    # Constants
    epsilon = math.radians(23.44)  # Obliquity of the ecliptic in radians

    # Calculate the number of days since the epoch J2000.0
    n = JD - 2451545.0

    # Calculate the mean longitude (L) and mean anomaly (g) of the Sun
    L = (280.46 + 0.9856474 * n) % 360  # Mean longitude in degrees
    g = math.radians((357.528 + 0.9856003 * n) % 360)  # Mean anomaly in radians

    # Calculate the ecliptic longitude (lambda) of the Sun
    lambda_sun = (L + 1.915 * math.sin(g) + 0.020 * math.sin(2 * g)) % 360

    # Convert lambda to radians for further calculations
    lambda_sun_rad = math.radians(lambda_sun)

    # Calculate the Right Ascension (alpha) in radians
    alpha_rad = math.atan2(math.cos(epsilon) * math.sin(lambda_sun_rad), math.cos(lambda_sun_rad))

    # Convert Right Ascension from radians to degrees
    alpha = math.degrees(alpha_rad)

    # Normalize alpha to be between 0 and 360 degrees
    alpha = alpha % 360

    # Convert Right Ascension from degrees to hours (since 360 degrees = 24 hours)
    RA_hours = alpha / 15.0

    return L, alpha, RA_hours  # Return mean longitude L in degrees, right ascension alpha in degrees, RA in hours

def calculate_eot(L, alpha):
    # Equation of Time in degrees
    EoT_deg = L - alpha

    # Ensure EoT is between -180 and +180 degrees
    if EoT_deg > 180:
        EoT_deg -= 360
    elif EoT_deg < -180:
        EoT_deg += 360

    # Convert EoT from degrees to minutes (since 1 degree = 4 minutes)
    EoT_min = EoT_deg * 4

    return EoT_min  # EoT in minutes

def calculate_solar_noon(longitude, EoT_min, timezone_offset):
    # longitude in degrees (positive east, negative west)
    # EoT_min in minutes
    # timezone_offset in hours (e.g., UTC+4 => timezone_offset = 4)
    SolarNoon = (12.0 - (longitude / 15.0) - (EoT_min / 60.0) + timezone_offset) % 24
    return SolarNoon  # in hours

def calculate_hour_angle(h, phi, delta):
    # h: solar altitude angle in degrees (negative if below the horizon)
    # phi: observer's latitude in degrees
    # delta: solar declination in degrees

    # Convert degrees to radians
    h_rad = math.radians(h)
    phi_rad = math.radians(phi)
    delta_rad = math.radians(delta)

    # Calculate cos(omega)
    cos_omega = (math.sin(h_rad) - math.sin(phi_rad) * math.sin(delta_rad)) / (math.cos(phi_rad) * math.cos(delta_rad))

    # Ensure cos_omega is within [-1, 1]
    cos_omega = min(1, max(-1, cos_omega))

    # Calculate omega in radians
    omega_rad = math.acos(cos_omega)

    # Convert omega to degrees
    omega = math.degrees(omega_rad)

    return omega  # Omega is between 0 and 180 degrees

def calculate_prayer_time(h, phi, delta, SolarNoon, before_noon=True):
    # h: solar altitude angle in degrees
    # phi: observer's latitude in degrees
    # delta: solar declination in degrees
    # SolarNoon: in hours

    # Calculate the Hour Angle omega
    omega = calculate_hour_angle(h, phi, delta)

    # Time difference in hours
    time_diff = omega / 15.0  # Earth rotates 15 degrees per hour

    # Determine if the event is before or after solar noon
    if before_noon:
        prayer_time = SolarNoon - time_diff
    else:
        prayer_time = SolarNoon + time_diff

    # Ensure prayer_time is between 0 and 24 hours
    prayer_time = prayer_time % 24

    return prayer_time

def calculate_prayer_times(Y, M, D, latitude, longitude, timezone_offset, asr_method='Standard', convention='MWL', is_ramadan=False):
    # Y, M, D: date
    # latitude, longitude: observer's location in degrees
    # timezone_offset: in hours (e.g., UTC+4 => timezone_offset = 4)
    # asr_method: 'Standard' for shadow ratio 1, 'Hanafi' for shadow ratio 2
    # convention: calculation convention for Fajr and Isha
    # is_ramadan: boolean indicating if it's Ramadan (affects Umm al-Qura calculation)

    hour, minute, second = 0, 0, 0  # Set time to midnight

    jd = calculate_julian_day(Y, M, D, hour, minute, second)

    declination = solar_declination(jd)

    L, alpha, RA_hours = calculate_right_ascension(jd)

    EoT_min = calculate_eot(L, alpha)

    SolarNoon = calculate_solar_noon(longitude, EoT_min, timezone_offset)

    # Define prayer angles based on the selected convention
    conventions = {
        'MWL': {'Fajr': 18.0, 'Isha': 17.0},  # Muslim World League
        'ISNA': {'Fajr': 15.0, 'Isha': 15.0},  # Islamic Society of North America
        'Egypt': {'Fajr': 19.5, 'Isha': 17.5},  # Egyptian General Authority of Survey
        'Makkah': {'Fajr': 18.5, 'Isha': None},  # Umm al-Qura University, Makkah
        'Karachi': {'Fajr': 18.0, 'Isha': 18.0},  # University of Islamic Sciences, Karachi
        'Tehran': {'Fajr': 17.7, 'Isha': 14.0},  # Institute of Geophysics, University of Tehran
    }

    if convention not in conventions:
        messagebox.showerror("Invalid Convention", "Selected convention is not recognized.")
        return

    fajr_angle = conventions[convention]['Fajr']
    isha_angle = conventions[convention]['Isha']

    # Define prayer angles in degrees
    prayer_angles = {
        'Fajr': -fajr_angle,      # Angle for Fajr
        'Sunrise': -0.833,  # Angle for Sunrise (includes atmospheric refraction)
        'Dhuhr': 0.0,       # Dhuhr at solar noon
        'Asr': None,        # Asr requires special calculation
        'Maghrib': -0.833,  # Angle for Maghrib (Sunset)
        'Isha': -isha_angle if isha_angle else None  # Angle for Isha or None if fixed time
    }

    prayer_times = {}

    # Fajr
    h = prayer_angles['Fajr']
    prayer_time = calculate_prayer_time(h, latitude, declination, SolarNoon, before_noon=True)
    prayer_times['Fajr'] = prayer_time

    # Sunrise
    h = prayer_angles['Sunrise']
    prayer_time = calculate_prayer_time(h, latitude, declination, SolarNoon, before_noon=True)
    prayer_times['Sunrise'] = prayer_time

    # Dhuhr
    prayer_times['Dhuhr'] = SolarNoon

    # Asr
    # Asr calculation depends on the shadow ratio
    if asr_method == 'Hanafi':
        asr_shadow_ratio = 2
    else:
        asr_shadow_ratio = 1  # Standard (Shafi'i, Maliki, Hanbali)

    phi_rad = math.radians(latitude)
    delta_rad = math.radians(declination)

    # Compute the angle for Asr
    # h = arccot(asr_shadow_ratio + tan(|phi - delta|))
    # arccot(x) = arctan(1 / x)
    denominator = asr_shadow_ratio + math.tan(abs(phi_rad - delta_rad))
    h_asr_rad = math.atan(1 / denominator)
    h_asr = math.degrees(h_asr_rad)

    # Now calculate the prayer time for Asr
    prayer_time = calculate_prayer_time(h_asr, latitude, declination, SolarNoon, before_noon=False)
    prayer_times['Asr'] = prayer_time

    # Maghrib
    h = prayer_angles['Maghrib']
    prayer_time = calculate_prayer_time(h, latitude, declination, SolarNoon, before_noon=False)
    prayer_times['Maghrib'] = prayer_time

    # Isha
    if convention == 'Makkah':
        # Umm al-Qura University, Makkah
        if is_ramadan:
            # 120 minutes after Maghrib during Ramadan
            isha_time = prayer_times['Maghrib'] + 2.0  # 2 hours
        else:
            # 90 minutes after Maghrib
            isha_time = prayer_times['Maghrib'] + 1.5  # 1.5 hours
        isha_time = isha_time % 24
        prayer_times['Isha'] = isha_time
    else:
        h = prayer_angles['Isha']
        prayer_time = calculate_prayer_time(h, latitude, declination, SolarNoon, before_noon=False)
        prayer_times['Isha'] = prayer_time

    return prayer_times

# Function to format time in HH:MM format
def format_time(time_in_hours):
    hours = int(time_in_hours)
    minutes = int(round((time_in_hours - hours) * 60))
    if minutes == 60:
        hours += 1
        minutes = 0
    return f"{hours:02d}:{minutes:02d}"

# Function to get user's IP-based location
def get_location():
    try:
        response = requests.get('https://ipapi.co/json/')
        data = response.json()
        latitude = float(data['latitude'])
        longitude = float(data['longitude'])
        city = data['city']
        return latitude, longitude, city
    except:
        return None, None, None

# Function to get timezone from latitude and longitude
def get_timezone(latitude, longitude):
    tf = TimezoneFinder()
    timezone_str = tf.timezone_at(lng=longitude, lat=latitude)
    return timezone_str

# Function to check if a date falls within Ramadan
def is_ramadan(date_obj):
    # For simplicity, define approximate Gregorian dates for Ramadan
    # Ramadan dates vary each year; for accurate calculation, you may need to use an Islamic calendar library
    # Here, we'll define Ramadan 2024 dates as an example
    ramadan_start = datetime(2024, 3, 11)
    ramadan_end = datetime(2024, 4, 9)
    return ramadan_start <= date_obj <= ramadan_end

# GUI Application
def calculate_and_display_prayer_times():
    # Get inputs from user
    date_str = date_entry.get()
    try:
        date_obj = datetime.strptime(date_str, '%Y-%m-%d')
    except ValueError:
        messagebox.showerror("Invalid Date", "Please enter the date in YYYY-MM-DD format.")
        return

    Y, M, D = date_obj.year, date_obj.month, date_obj.day

    try:
        latitude = float(lat_entry.get())
        longitude = float(lon_entry.get())
    except ValueError:
        messagebox.showerror("Invalid Coordinates", "Please enter valid numbers for latitude and longitude.")
        return

    asr_method = asr_var.get()

    # Get timezone offset
    timezone_str = get_timezone(latitude, longitude)
    if timezone_str:
        timezone = pytz.timezone(timezone_str)
        timezone_offset = timezone.utcoffset(datetime.now()).total_seconds() / 3600.0
    else:
        try:
            timezone_offset = float(timezone_entry.get())
        except ValueError:
            messagebox.showerror("Invalid Timezone", "Please enter a valid number for timezone offset.")
            return

    convention = convention_var.get()

    # Check if it's Ramadan for Umm al-Qura convention
    ramadan = False
    if convention == 'Makkah':
        ramadan = is_ramadan(date_obj)

    # Calculate prayer times
    prayer_times = calculate_prayer_times(Y, M, D, latitude, longitude, timezone_offset, asr_method, convention, ramadan)

    # Display prayer times
    result = f"Prayer times for {date_str} at latitude {latitude}°, longitude {longitude}°:\n"
    for prayer, time in prayer_times.items():
        result += f"{prayer}: {format_time(time)}\n"
    messagebox.showinfo("Prayer Times", result)

def auto_fill_location():
    latitude, longitude, city = get_location()
    if latitude and longitude:
        lat_entry.delete(0, tk.END)
        lat_entry.insert(0, f"{latitude:.6f}")
        lon_entry.delete(0, tk.END)
        lon_entry.insert(0, f"{longitude:.6f}")
        location_label.config(text=f"Detected Location: {city}")
        # Get timezone
        timezone_str = get_timezone(latitude, longitude)
        if timezone_str:
            timezone_entry.delete(0, tk.END)
            timezone = pytz.timezone(timezone_str)
            timezone_offset = timezone.utcoffset(datetime.now()).total_seconds() / 3600.0
            timezone_entry.insert(0, f"{timezone_offset}")
    else:
        messagebox.showwarning("Location Error", "Could not detect location automatically.")

# Main window
root = tk.Tk()
root.title("Prayer Times Calculator")

# Date input
tk.Label(root, text="Date (YYYY-MM-DD):").grid(row=0, column=0, sticky=tk.W)
date_entry = tk.Entry(root)
date_entry.grid(row=0, column=1)
date_entry.insert(0, datetime.now().strftime('%Y-%m-%d'))

# Latitude input
tk.Label(root, text="Latitude:").grid(row=1, column=0, sticky=tk.W)
lat_entry = tk.Entry(root)
lat_entry.grid(row=1, column=1)

# Longitude input
tk.Label(root, text="Longitude:").grid(row=2, column=0, sticky=tk.W)
lon_entry = tk.Entry(root)
lon_entry.grid(row=2, column=1)

# Timezone offset input
tk.Label(root, text="Timezone Offset (hours):").grid(row=3, column=0, sticky=tk.W)
timezone_entry = tk.Entry(root)
timezone_entry.grid(row=3, column=1)

# Asr method selection
tk.Label(root, text="Asr Calculation Method:").grid(row=4, column=0, sticky=tk.W)
asr_var = tk.StringVar(value='Standard')
tk.Radiobutton(root, text="Standard (Shadow Ratio 1)", variable=asr_var, value='Standard').grid(row=4, column=1, sticky=tk.W)
tk.Radiobutton(root, text="Hanafi (Shadow Ratio 2)", variable=asr_var, value='Hanafi').grid(row=5, column=1, sticky=tk.W)

# Convention selection
tk.Label(root, text="Calculation Convention:").grid(row=6, column=0, sticky=tk.W)
convention_var = tk.StringVar(value='MWL')
conventions_list = [
    ('Muslim World League', 'MWL'),
    ('Islamic Society of North America (ISNA)', 'ISNA'),
    ('Egyptian General Authority of Survey', 'Egypt'),
    ('Umm al-Qura University, Makkah', 'Makkah'),
    ('University of Islamic Sciences, Karachi', 'Karachi'),
    ('Institute of Geophysics, Tehran', 'Tehran'),
]
row = 6
for text, value in conventions_list:
    tk.Radiobutton(root, text=text, variable=convention_var, value=value).grid(row=row, column=1, sticky=tk.W)
    row += 1

# Auto-fill location button
auto_location_button = tk.Button(root, text="Auto-detect Location", command=auto_fill_location)
auto_location_button.grid(row=row, column=0, columnspan=2, pady=5)
row += 1

# Detected location label
location_label = tk.Label(root, text="")
location_label.grid(row=row, column=0, columnspan=2)
row += 1

# Calculate button
calculate_button = tk.Button(root, text="Calculate Prayer Times", command=calculate_and_display_prayer_times)
calculate_button.grid(row=row, column=0, columnspan=2, pady=10)

root.mainloop()
