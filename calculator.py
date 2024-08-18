import math

def calculate_julian_day(Y, M, D, hour=0, minute=0, second=0):
    if M <= 2:
        Y -= 1
        M += 12

    A = math.floor(Y / 100)
    B = 2 - A + math.floor(A / 4)

    # Calculate the Julian Day Number (JDN)
    JD = (math.floor(365.25 * (Y + 4716)) +
          math.floor(30.6001 * (M + 1)) +
          D + B - 1524.5)

    # Convert time to a fraction of a day and add to JD
    day_fraction = (hour + minute / 60 + second / 3600) / 24.0
    JD += day_fraction

    return JD

# Example usage:
Y, M, D = 2024, 8, 18
hour, minute, second = 15, 30, 0  # 3:30 PM
jd = calculate_julian_day(Y, M, D, hour, minute, second)
print(f"The Julian Day for {Y}-{M}-{D} {hour}:{minute}:{second} is {jd}")

def solar_declination(JD):
    # JD for January 1st of the year corresponding to the input JD
    JD_2000 = 2451545.0  # Julian Day for January 1, 2000, at 12:00 PM (noon)
    days_since_2000 = JD - JD_2000

    # Number of days in the year
    N = days_since_2000 % 365.24
    
    # Calculate solar declination in radians
    declination_radians = math.asin(math.sin(math.radians(-23.44)) * math.cos(math.radians((360 / 365.24) * (N + 10))))
    
    # Convert to degrees
    declination_degrees = math.degrees(declination_radians)
    
    return declination_degrees

#usage:
declination = solar_declination(jd)
print(f"Solar declination for JD {jd} is {declination:.2f} degrees")