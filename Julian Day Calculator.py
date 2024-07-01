import datetime

# Input year, month, and day as ints
Y = int(input("Input the year: "))
M = int(input("Input the month: "))
D = int(input("Input the day: "))

# Input universal time
time = input("Input the universal time (HH:MM:SS): ")

# Parse the time string
time_obj = datetime.datetime.strptime(time, "%H:%M:%S")

# Extract hours, minutes, and seconds
hours, minutes, seconds = time_obj.hour, time_obj.minute, time_obj.second

# Calculate UT in decimal hours
UT = hours + (minutes/60) + (seconds/3600)

# Calculate Julian Date
if M <= 2:
    Y -= 1
    M += 12

A = int(Y / 100)
B = 2 - A + int(A / 4)

JD = int(365.25 * (Y + 4716)) + int(30.6001 * (M + 1)) + D + B - 1524.5 + (UT / 24)

print(f"The Julian Date is: {JD}")