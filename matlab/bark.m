function z = bark(freq)

z = 13*atan(0.76*(freq/1000)) + 3.5*atan((freq/7500).^2);