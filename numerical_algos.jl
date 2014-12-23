# Copyright (c) 2015, Valentino Constantinou <vc1492a@gmail.com>
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

#distance between two points. Can easily extend to R^n.
function point_distance_R2(x1, y1, x2, y2)
    distance = sqrt(((x2 - x1) ^ 2) + ((y2 - y1) ^ 2))
    return distance
end

function point_distance_R3(x1, y1, z1, x2, y2, z2)
    distance = sqrt(((x2 - x1) ^ 2) + ((y2 - y1) ^ 2) + ((z2 - z1) ^ 2))
    return distance
end

point_distance_R2(0., 1., 1., 0.) #Pythagorean Theorem tells us this is the correct result.
point_distance_R3(1., 2., 3., 3., 2., 1.) #correct according to Wolfram Alpha

#calculate the hypotenuse of a triangle
function hypot(x,y)
  x = abs(x)
  y = abs(y)
  if x > y
    r = y/x
    return x*sqrt(1+r*r)
  end
  if y == 0
    return zero(x)
  end
  r = x/y
  return y*sqrt(1+r*r)
end

hypot(3,4)

#area of a triangle using points
function tri_point_area(x1, x2, x3, y1, y2, y3)
    A = [x1 y1 1; x2 y2 1; x3 y3 1]
    area = 0.5 * abs(det(A))
    return area
end

tri_point_area(-1, 0, 1, 0, 1, 0)

#real and absolute error
function absError(real, approx)
    abserror = (approx - real)
    return abserror
end

function realError(real, approx)
    realerror = approx / real
    return realerror
end

#trapezoid method
function trapezoid_uniform(f, a, b, n)
    h = (b - a) / n
    s = .5 * (f(a) + f(b))
    for i =1:(n - 1)
        x = a + i * h
        s = s + f(x)
    end
    return s * h
end

function f1(x)
    fn = sin(x)
    return fn
end

trapezoid_uniform(f1, 0, pi, 100) #returns 1.9998355038874436. 

#trapazoid method, while generally affective, does not estimate functions like f(x)=e^x well. This is due to the
#exponentially increasing nature of the function. The higher the magnitude of the slope, the more error we can
#expect from using the trapezoid method.

#simpson's rule
function simpson(f, a, b, level, level_max, precision)
    level += 1
    h = b - a
    c = (a + b) / 2
    one_s = h * (f(a) + 4. * f(c) + f(b)) / 6.
    d = (a + c) / 2
    e = (c + b) / 2
    two_s = h * (f(a) + 4. * f(d) + 2. * f(c) + 4. * f(e) + f(b)) / 12.
    if level >= level_max
        simpson_result = two_s
        print("Max level reached")
        return simpson_result
    else
        if abs(two_s - one_s) < 15. * precision
            return two_s + (two_s - one_s) / 15.
        else
            left_s = simpson(f, a, c, level, level_max, precision / 2.)
            right_s = simpson(f, c, b, level, level_max, precision / 2.)
            print(left_s)
            return left_s + right_s
        end
    end
end

function simptest(x)
    fn = x^2
    return fn
end

simpson(simptest, -2, 2, 0.000000000000001, 1, 10) #returns 5.333333333333.

#open Newton-Cotes rules
function nc_midpoint(f, a, b, n)
    h = 0.5 * (b - a) #midpoint
    x1 = a + h #we need to move one step away from a. 
    sum = 2. * h * f(x1) #formula
    return sum
end

function nc_two_point(f, a, b, n)
    h = 0.3333333 * (b - a) #divide interval by 3
    x1 = a + h
    x2 = a + 2. * h #we need to move twice the step size away from a.
    sum = 1.5 * h * (f(x1) + f(x2)) #formula
    return sum
end

function nc_three_point(f, a, b, n)
    h = 0.25 * (b - a) #divide interval by 4
    x1 = a + h #we need to move one step away from a. 
    x2 = a + 2. * h #we need to move twice the step size away from a.
    x3 = a + 3. * h #we need to move three times the step size away from a.
    sum = 1.33333333 * h * (2. * f(x1) - f(x2) + 2. * f(x3)) #formula
    return sum
end

function nc_four_point(f, a, b, n)
    h = 0.2 * (b - a) #divide interval by 5
    x1 = a + h #we need to move one step away from a. 
    x2 = a + 2. * h #we need to move twice the step size away from a.
    x3 = a + 3. * h #we need to move three times the step size away from a.
    x4 = a + 4. * h #we need to move four times the step size away from a.
    sum = 0.208333333 * h* (11. * f(x1) + f(x2) + f(x3) + 11. * f(x4)) #formula
    return sum
end

function nc_five_point(f, a, b, n)
    h = 0.1666666667 * (b - a) #divide interval by 6
    x1 = a + h #we need to move one step away from a. 
    x2 = a + 2. * h #we need to move twice the step size away from a.
    x3 = a + 3. * h #we need to move three times the step size away from a.
    x4 = a + 4. * h #we need to move four times the step size away from a.
    x5 = a + 5. * h #we need to move five times the step size away from a.
    sum = 0.3 * h * (11. * f(x1) - 14. * f(x2) + 26. * f(x3) - 14. * f(x4) + 11. * f(x5)) #formula
    return sum
end


#Problem 5.3.7

#let's test the open Newton-Cotes rules. 

function f(x) #to calculate the above function, we need an open interval (-1,1). There are asymptotes -1 and 1.
    fn = 1 /((1 - (x^2)) ^ (0.5))
    return fn
end

nc_midpoint(f, -1, 1, 10) #returns 2.0
nc_two_point(f, -1, 1, 10) #returns 2.121320104911125
nc_three_point(f, -1, 1, 10) #returns 2.412534762980001
nc_four_point(f, -1, 1, 10) #returns 2.4617701170877777
nc_five_point(f, -1, 1, 10) #returns 2.581761250230586

#The true result is 3.14159, or pi. You can see here that even with the five-point rule, we are still far away from
#the true result. This is due to the nature of the function, which as vertical asymptotes at x=-1 and x=1.
#Theoretically, adding additional points should improve the precision of the result, but this is impractical
#given other methods of integration, such as Gaussian quadrature or Richardson extrapolation. 

#Gaussian Quadrature
#Two point Gaussian Quadrature
function two_pt_gauss(f, a, b, n)
    h = (b - a) / n
    #print h
    sum = 0
    for i = 1:n-1
        x0 = a + (i * h) #starting at left end point, h represents step size.
        #print x0
        x1 = x0 + (0.5 * h) * (1 - sqrt(1 / 3))
        #print x1
        x2 = x0 + (0.5 * h) * (1 + sqrt(1 / 3))
        #print x2
        sum += ((f(x1) + f(x2))) #weights are 1.
        #print sum
    end
    sum *= (0.5 * h)
    return sum
end

#Three point Gaussian Quadrature
function three_pt_gauss(f, a, b, n)
    h = (b - a) / n
    #print h
    sum = 0
    for i = 1:n-1
        x0 = a + (i * h) #starting at left end point, h represents step size.
        #print x0
        x1 = x0 + (0.5 * h) * (1 - sqrt(3 / 5))
        #print x1
        x2 = x0 + (0.5 * h)
        #print x2
        x3 = x0 + (0.5 * h) * (1 + sqrt(3 / 5))
        #print x3
        sum += (((5 / 9) * f(x1)) + ((8 / 9) * f(x2)) + ((5 / 9) * f(x3))) #weights are 5/9, 8/9, and 5/9.
        #print sum
    end
    sum *= (0.5 * h)
    return sum
end


function f(x)
    fn = x^5
    return fn
end

two_pt_gauss(f,0.,1.,10) #result is 0.16666512500000005
three_pt_gauss(f,0.,1.,10) #result is 0.16666650000000005

#You can see here that using Gaussian Quadrature feels like cheating, almost. With so little computation,
#we can approximate the area under the curve (i.e. integrate) with a, compared to other methods, very high
#amount of accuracy. Using n=10 and three_pt_gauss, we obtain the true result.

