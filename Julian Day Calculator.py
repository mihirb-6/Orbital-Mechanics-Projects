# %%
import datetime

print("Input the year:")
Y = float(input())
print("Input the month:")
M = float(input())
print("Input the day:")
D = float(input())
print("Input the universal time:")
time = str(input())

time_obj = datetime.datetime.strptime(time, "%H:%M:%S")

hour = time_obj.hour
minute = time_obj.minute
seconds = time_obj.second
# %%
UT = hour + (minute/60) + (seconds/3600)

J0 = (367*Y) - int((7*(Y+int((M+9)/12)))/4) + int((275*M)/9)
+ D + 17121013.5

print(UT)
print(J0)
JD = J0 + (UT / 24)

print(f"The Julian Day Number is: {JD}")
