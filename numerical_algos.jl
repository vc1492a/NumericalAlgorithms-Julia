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

#Function for creating a Hilbert Matrix
function createHilbert(n) #define a function for creating Hilbert matrices.
    H = zeros(n, n) #declare an empty array of nxn dimension.
    for i = 1:n
        for j = 1:n
            H[i,j] = 1 / (i + j - 1) #formula for constructing matrix values.
        end
    end
    return H #return the matrix A
end

createHilbert(5)

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

#Newton's Method
function newton(f,f_prime,x,tolerance,precision,display,steps)
    fx = f(x)
    for i = 1:steps
        if abs(f_prime(x)) < precision
            print("small derivative")
            break
        end
        d = fx / f_prime(x)
        x = x - d
        fx = f(x)
        if display == 1
                print([i, x, fx])
        end
        if abs(d) < tolerance
            print("convergence")
            break
        end
    end
end

function newtonaccel(f,f_prime,x,tolerance,precision,display,steps)
    fx = f(x)
    for i = 1:steps
        if abs(f_prime(x)) < precision
            print("small derivative")
            break
        end
        d = fx / f_prime(x)
        x = x - d
        fx = f(x)
        if display == 1
                print([i, x, fx])
        end
        if abs(d) < tolerance
            print("convergence")
            break
        end
    end
end

function newtonmod(f,f_prime,x,tolerance,precision,display,steps)
    fx = f(x)
    for i = 1:steps
        if abs(f_prime(x)) < precision
            print("small derivative")
            break
        end
        d = fx / f_prime(x)
        x = x - d
        if abs(f(x - d)) >= abs(f(x))
            d = 0.5 * d
        else
            d = d
        end
        fx = f(x)
        if display == 1
                print([i, x, fx])
        end
        if abs(d) < tolerance
            print("convergence")
            break
        end
    end
end

function testfn(x)
    fn = sin(x)
    return fn
end

function testfnprime(x)
    fn = cos(x)
    return fn
end

newton(testfn, testfnprime, 1.2, 0.00000001, 0.00000001, 1, 25)
newtonaccel(testfn, testfnprime, 1.2, 0.00000001, 0.00000001, 1, 25)
newtonmod(testfn, testfnprime, 1.2, 0.00000001, 0.00000001, 1, 25)

#secant method
function secant(f, a, b, precision, steps)
    fa = f(a)
    fb = f(b)
    if abs(fa) > abs(fb)
        temp = a
        a = b
        b = temp
        temp = fa
        fa = fb
        fb = temp
    end
    print([0,a,fa])
    print([1,b,fb])
    for n = 2:steps
        if abs(fa) > abs(fb)
            temp = a
            a = b
            b = temp
            temp = fa
            fa = fb
            fb = temp
        end
        d = (b - a) / (fb - fa)
        b = a
        fb = fa
        d = d * fa
        if abs(d) < precision
            print("Convergence")
            print([n, a, fa])
        end
        a = a - d
        fa = f(a)
    end
end

function secantfn(x)
    fn = (e^x)-3*(x^2)
    return fn
end

secant(secantfn, -0.5, 2. ,1. ^ -10, 25)
#found root to be -0.45896257524. There are other roots at 0.91 and 3.733

#simpson's Method
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

#Chebyshev Nodes
function chebyshevnodes(func,a,b,n) #number of nodes defined by user as n, a=left most point and b=right most point.
    x = Float64[] #empty matrix to store values.
    y = Float64[]
    for i = 1 #really only need to check this once.
        if 0<=i<=n #diagnostic to make sure we are operating within the bounds of the function. Not really needed.
            print("i is less than or equal to n")
        else
            print("i is NOT less than or equal to n")
            break
        end
    end
    for i = 1:n+1
        push!(x, cos((((2.0*float(i)+1)/(2.0*float(n)+2))*pi))) #formula given in book.
        push!(y, func(x[i]))
    end
    return [x,y] #returns calculated chebyshev nodes to matrix for use.
end

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

#Runge-Kutta methods
# The Runge-Kutta method imitate the Taylor series method without requiring analytic differentiation of the
# original differential equation. Therefore, the below algorithms will accept any function, interval, and given point
# and will represent the solution to an ODE locally at the given point. This method works provided we now the
# value exactly at some arbitrary t. The downside of these methods is that they have to evaluate the function f several times,
# which can be very time consuming depending on the function.

#f denotes the function, x denotes the initial point, a and b denote the beginning and end of the interval, n denotes number of iterations.

# Runge-Kutta of order 2.
function runge_kutta_2(f, x, a, b, n)
    h = (b - a) / n  # step size
    t = a  # sets the initial t as the left-most point of the interval, a.
    for j = 1:n
        k1 = h * f(t, x)
        k2 = h * f(t + h, x + k1)
        x += 0.5 * (k1 + k2)
        t = a + (j * h)
        print([j, t, x])
    end
end

# Runge-Kutta method of order 4.
function runge_kutta_4(f, x, a, b, n)
    h = (b - a) / n  # step size
    t = a  # sets the initial t as the left-most point of the interval, a.
    for j = 1:n
        k1 = h * f(t, x)
        k2 = h * f(t + 0.5 * h, x + 0.5 * k1)
        k3 = h * f(t + 0.5 * h, x + 0.5 * k2)
        k4 = h * f(t + h, x + k3)
        x += (1. / 6.) * (k1 + 2. * k2 + 2. * k3 + k4)
        t = a + (j * h)
        print([j, t, x])
    end
end

#For the Runge-Method of order 5 given below,
#the difference between the values of x(t+h) obtained from the 4th and 5th order procedures is an estimate of the
#local truncation error in the 4th order procedure. Therefore, evaluations give a 5th order approximation,
#together with an error estimate. This method is called the Runge-Kutta-Fehlberg method

#adaptive (enter 1 for yes, 0 for no) is a toggle to store error values needed in the adaptive Runge-Kutta method.

# Runge-Kutta method of order 5.
function runge_kutta_fehlberg(f, x, a, b, n, adaptive)
    h = (b - a) / n  # step size
    t = a  # sets the initial t as the left-most point of the interval, a.
    if adaptive == 1
        erstore = Float64[] #need a place to store error values for adaptive Runge-Kutta.
    end
    #note that a2 == b2 == 0.
    c20 = 0.25
    c21 = 0.25
    c30 = 0.375
    c31 = 0.09375
    c32 = 0.28125
    c40 = 12. / 13.
    c41 = 1932. / 2197.
    c42 = 7200. / 2197.
    c43 = 7296. / 2197.
    c51 = 439. / 216.
    c52 = -8.
    c53 = 3680. / 513.
    c54 = -845. / 4104.
    c60 = 0.5
    c61 = -8. / 27.
    c62 = 2.
    c63 = -3544. / 2565.
    c64 = 1859. / 4104.
    c65 = -0.275
    a1 = 25. / 216.
    a3 = 1408. / 2565.
    a4 = 2197. / 4104.
    a5 = -0.2
    b1 = 16. / 135.
    b3 = 6656. / 12825.
    b4 = 28561. / 56430.
    b5 = -0.18
    b6 = 2. / 55.
    for j = 1:n
        k1 = h * f(t, x)
        k2 = h * f(t + (c20 * h), x + (c21 * k1))
        k3 = h * f(t + (c30 * h), x + (c31 * k1) + (c32 * k2))
        k4 = h * f(t + (c40 * h), x + (c41 * k1) + (c42 * k2) + (c43 * k3))
        k5 = h * f(t + h, x + (c51 * k1) + (c52 * k2) + (c53 * k3) + (c54 * k4))
        k6 = h * f(t + (c60 * h), x + (c61 * k1) + (c62 * k2) + (c63 * k3) + (c64 * k4) + (c65 * k5))
        x4 = x + ((a1 * k1) + (a3 * k3) + (a4 * k4) + (a5 * k5))
        x +=((b1 * k1) + (b3 * k3) + (b4 * k4) + (b5 * k5) + (b6 * k6))
        t = a + (j * h)
        er = abs(x - x4)
        if adaptive == 1
            push!(erstore, er)
        else
            print([j, t, x, er])
        end
    end
    if adaptive == 1
        return erstore
    end
end

#The Bogacki Shampine method is a Runge Kutta method of order three with four stages with the First
#Same As Last (FSAL) property, so that it uses approximately three function evaluations per step.
#It has an embedded second-order method which can be used to implement adaptive step size.

function runge_kutta_bogacki_shampine(f, x, a, b, n, adaptive)
    h = (b - a) / n  # step size
    t = a  # sets the initial t as the left-most point of the interval, a.
    if adaptive == 1
        erstore = Float64[] #need a place to store error values for adaptive Runge-Kutta.
    end
    #note that b4 == c31 == 0
    c20 = 0.5
    c21 = 0.5
    c30 = 0.75
    c32 = 0.75
    c40 = 1.
    c41 = 2. / 9.
    c42 = 1. / 3.
    c43 = 4. / 9.
    a1 = 2. / 9.
    a2 = 1. / 3.
    a3 = 4. / 9.
    b1 = 7. / 24.
    b2 = 1. / 4.
    b3 = 1. / 3.
    b4 = 1. / 8.
    for j in 1:n
        k1 = h * f(t, x)
        k2 = h * f(t + (c20 * h), x + (c21 * k1))
        k3 = h * f(t + (c30 * h), x + (c32 * k2))
        k4 = h * f(t + (c40 * h), x + (c41 * k1) + (c42 * k2) + (c43 * k3))
        x3 = x + ((a1 * k1) + (a2 * k2) + (a3 * k3))
        x +=((b1 * k1) + (b2 * k2) + (b3 * k3) + (b4 * k4))
        t = a + (j * h)
        er = abs(x - x3)
        if adaptive == 1
            push!(erstore, er)
        else
            print([j, t, x, er])
        end
    end
    if adaptive == 1
        return erstore
    end
end

#The Cash-Karp Runge-Kutta method uses six function evaluations to calculate fourth and fifth order accurate solutions.
#The difference between these solutions is then taken to be the error of the fourth order solution.

function runge_kutta_cash_karp(f, x, a, b, n, adaptive)
    h = (b - a) / n  # step size
    t = a  # sets the initial t as the left-most point of the interval, a.
    if adaptive == 1
        erstore = Float64[] #need a place to store error values for adaptive Runge-Kutta.
    end
    #note that a2 == a5 == b2 == 0.
    c20 = 1. / 5.
    c21 = 1. / 5.
    c30 = 3./ 10.
    c31 = 3./ 40.
    c32 = 9./ 40.
    c40 = 3. / 5.
    c41 = 3. / 10.
    c42 = -9. / 10.
    c43 = 6. / 5.
    c50 = 1.
    c51 = -11. / 54.
    c52 = 5. / 2.
    c53 = -70. / 27.
    c54 = 35. / 27.
    c60 = 7. / 8.
    c61 = 1631. / 55296.
    c62 = 175. / 512.
    c63 = 575. / 13824.
    c64 = 44275. / 110592.
    c65 = 253. / 2096.
    a1 = 37. / 378.
    a3 = 250. / 621.
    a4 = 125. / 594.
    a6 = 512. / 1771.
    b1 = 2825. / 27648.
    b3 = 18575. / 48384.
    b4 = 13525. / 55296.
    b5 = 277. / 14336.
    b6 = 1. / 4.
    for j = 1:n
        k1 = h * f(t, x)
        k2 = h * f(t + (c20 * h), x + (c21 * k1))
        k3 = h * f(t + (c30 * h), x + (c31 * k1) + (c32 * k2))
        k4 = h * f(t + (c40 * h), x + (c41 * k1) + (c42 * k2) + (c43 * k3))
        k5 = h * f(t + (c50 * h), x + (c51 * k1) + (c52 * k2) + (c53 * k3) + (c54 * k4))
        k6 = h * f(t + (c60 * h), x + (c61 * k1) + (c62 * k2) + (c63 * k3) + (c64 * k4) + (c65 * k5))
        x4 = x + ((a1 * k1) + (a3 * k3) + (a4 * k4) + (a6 * k6))
        x += (b1 * k1) + (b3 * k3) + (b4 * k4) + (b5 * k5) + (b6 * k6)
        t = a + (j * h)
        er = abs(x - x4)
        if adaptive == 1
            push!(erstore, er)
        else
            print([j, t, x, er])
        end
    end
    if adaptive == 1
        return erstore
    end
end

#The Dormand Prince method is another type of ODE solver in the Runge-Kutta family.
#It uses six evaluations to calculate fourth and fifth order accurate solutions.
#The Dormand Pricnce has seven stages, but it uses only six function evaluations per step because
#it has the FSAL (First Same As Last) property: the last stage is evaluated at the same point as the first stage
#of the next step. Dormand and Prince chose the coefficients of their method to minimize the error of the fifth-order
#solution. This is the main difference with the Fehlberg method, which was constructed so that the fourth-order
#solution has small error. For this reason, the Dormand Prince method is more suitable when the higher-order
#solution is used to continue the integration.

function runge_kutta_dormand_prince(f, x, a, b, n, adaptive)
    h = (b - a) / n  # step size
    t = a  # sets the initial t as the left-most point of the interval, a.
    if adaptive == 1
        erstore = Float64[] #need a place to store error values for adaptive Runge-Kutta.
    end
    #note that a2 == a6 == b2 == b6 == c72 == 0.
    c20 = 1. / 5.
    c21 = 1. / 5.
    c30 = 3. / 10.
    c31 = 3. / 40.
    c32 = 9. / 40.
    c40 = 4. / 5.
    c41 = 44. / 45.
    c42 = -56. / 15.
    c43 = 32. / 9.
    c50 = 8. / 9.
    c51 = 19372. / 6561.
    c52 = -25360. / 2187.
    c53 = 64448. / 6561.
    c54 = -212. / 729.
    c60 = 1.
    c61 = 9017. / 3168.
    c62 = -355. / 33.
    c63 = 46732. / 5247.
    c64 = 49. / 176.
    c65 = -5103. / 18656.
    c70 = 1.
    c71 = 35. / 384.
    c73 = 500. / 1113.
    c74 = 125. / 132.
    c75 = -2187. / 6784.
    c76 = 11. / 84.
    a1 = 35. / 384.
    a3 = 500. / 1113.
    a4 = 125. / 192.
    a5 = -2187. / 6784.
    b1 = 5179. / 57600.
    b3 = 7571. / 16695.
    b4 = 393. / 640.
    b5 = -92097. / 339200.
    b6 = 187. / 2100.
    b7 = 1. / 40.
    for j = 1:n
        k1 = h * f(t, x)
        k2 = h * f(t + (c20 * h), x + (c21 * k1))
        k3 = h * f(t + (c30 * h), x + (c31 * k1) + (c32 * k2))
        k4 = h * f(t + (c40 * h), x + (c41 * k1) + (c42 * k2) + (c43 * k3))
        k5 = h * f(t + (c50 * h), x + (c51 * k1) + (c52 * k2) + (c53 * k3) + (c54 * k4))
        k6 = h * f(t + (c60 * h), x + (c61 * k1) + (c62 * k2) + (c63 * k3) + (c64 * k4) + (c65 * k5))
        k7 = h * f(t + (c70 * h), x + (c71 * k1) + (c73 * k3) + (c74 * k4) + (c75 * k5) + (c76 * k6))
        x5 = x + ((a1 * k1) + (a3 * k3) + (a4 * k4) + (a5 * k5))
        x += ((b1 * k1) + (b3 * k3) + (b4 * k4) + (b5 * k5) + (b6 * k6) + (b7 * k7))
        t = a + (j * h)
        er = abs(x - x5)
        if adaptive == 1
            push!(erstore, er)
        else
            print([j, t, x, er])
        end
    end
    if adaptive == 1
        return erstore
    end
end

#The error estimate 'er' calculated and stored in some of the above methods
#can tell us when to adjust the step size to control for the single-step error.
#This fact, together with the 5th order approximation, results in an adaptive procedure
#that is very accurate. emin and emax are the lower and upper bounds on the allowable error estimate.
#hmin and hmax are are bounds on the step size h.
#eflag is an error flag that returns 0 or 1 depending on if there is a successful march from a to b or if max. n is reached.
#method is a toggle to select between different Runge-Kutta methods, as follows.

#1: Runge-Kutta-Fehlberg
#2: Runge-Kutta Bogacki-Shampine
#3: Runge-Kutta Cash-Karp
#4: Runge-Kutta Dormand-Prince

function runge_kutta_adaptive(f, x, a, b, n, emin, emax, hmin, hmax, method)
    erstore = Float64[]
    if method == 1
        append!(erstore, runge_kutta_fehlberg(f, x, a, b, n, 1)) #store error values for use in adaptive method.
    end
    if method == 2
        append!(erstore, runge_kutta_bogacki_shampine(f, x, a, b, n, 1)) #store error values for use in adaptive method.
    end
    if method == 3
        append!(erstore, runge_kutta_cash_karp(f, x, a, b, n, 1)) #store error values for use in adaptive method.
    end
    if method == 4
        append!(erstore, runge_kutta_dormand_prince(f, x, a, b, n, 1)) #store error values for use in adaptive method.
    end
    h = (b - a) / n
    t = a
    eflag = 1
    k = 0
    sig = 0.5 * (10.0 ^ -5.)
    while k <= n
        #maybe add for loop here for erstore at the bottom of the function.
        k = k + 1
        if abs(h) < hmin
            h = copysign(hmin,(h)) #return h with sign(h)*hmin
        end
        if abs(h) > hmax
            h = copysign(hmax,(h)) #return h with sign(h)*hmax
        end
        d = abs(b - a)
        if d <= abs(h)
            eflag = 0
            if d <= sig * max(abs(b), abs(a))
                break
            end
            h = copysign(d,(h)) #return h with sign(h)*d
        end
        xsave = x
        tsave = t
        if method == 1
            runge_kutta_fehlberg(f, x, a, b, n, 0) #note: no need to return error values for printing results
        end
        if method == 2
            runge_kutta_bogacki_shampine(f, x, a, b, n, 0)
        end
        if method == 3
            runge_kutta_cash_karp(f, x, a, b, n, 0)
        end
        if method == 4
            runge_kutta_dormand_prince(f, x, a, b, n, 0)
        end
        if eflag == 0
            break
        end
        for i = 1:n
            if i in erstore < emin
                h = 2. * h
            end
            if i in erstore > emax
                h = 0.5 * h
                x = xsave
                t = tsave
                k = k - 1
            end
        end
    end
end


# test function given in book.
function f(t, x)
    fn = 2. + ((x - t - 1.) ^ 2.)
    return fn
end

#another test function
function g(t,x)
    fn = 3. + (5 * sin(t)) + (0.2 * x)
    return fn
end

#tests using f(t,x):
runge_kutta_2(f, 2., 1., 1.5625, 72) #returns 3.192942728232579
runge_kutta_4(f, 2., 1., 1,5625, 72) #returns 3.192937673837072
runge_kutta_fehlberg(f, 2., 1., 1.5625, 72, 0) #returns 3.2206213352344504
runge_kutta_bogacki_shampine(f, 2., 1., 1.5625, 72, 0) #returns 3.192940384731455
runge_kutta_cash_karp(f, 2., 1., 1.5625, 72, 0) #returns 3.193038424315954
runge_kutta_dormand_prince(f, 2., 1., 1.5625, 0) #returns 3.1929887485024313
runge_kutta_adaptive(f, 2., 1., 1.5625, 72, 10. ^ -8, 10. ^ -5, 10. ^ -6, 1.0, 4) #returns 3.1929887485024313

#tests using g(t,x):
runge_kutta_2(g, 0., 0., 10., 1000) #returns 135.91667314003618
runge_kutta_4(g, 0., 0., 10., 1000) #returns 135.91724460992168
runge_kutta_fehlberg(g, 0., 0., 10., 100000, 0) #returns 135.93915664626414
runge_kutta_bogacki_shampine(g, 0., 0., 10., 100000, 0) #returns 135.917244604469
runge_kutta_cash_karp(g, 0., 0., 10., 100000, 0) #returns 135.91734192973425
runge_kutta_dormand_prince(g, 0., 0., 10., 100000, 0) #returns 135.91729347356193
runge_kutta_adaptive(g, 0., 0., 10., 1000, 10. ^ -8, 10. ^ -5, 10. ^ -6, 1.0, 4) #returns 135.9221293554806

#Convert polar coordinates to cartesian coordinates
function polars_to_cartesian(r, d)
    theta = (d * pi) / 180.
    x = r*cos(theta)
    y = r*sin(theta)
    return x, y
end

polars_to_cartesian(3, 30)

#Jacobian Matrix evaluation
#given two (or three) functions with respect to x and y (and z), calculates the determinant of the Jacobian matrix at a given point.
#useful in characterizing local behavior of nonlinear system about an equilibrium point, i.e. linearization.

#fx,fy,gx,gy are functions f and g with respect to x and y
#xi,yi are points of which to evaluate the Jacobian

function jacobian_eval_R2(fx,fy,gx,gy,xi,yi)
    a = fx(xi,yi)
    b = fy(xi,yi)
    c = gx(xi,yi)
    d = gy(xi,yi)
    Jcomplete = [a b;c d]
    return det(Jcomplete)
end

function jacobian_eval_R3(fx,fy,fz,gx,gy,gz,hx,hy,hz,xi,yi,zi)
    a = fx(xi,yi,zi)
    b = fy(xi,yi,zi)
    c = fz(xi,yi,zi)
    d = gx(xi,yi,zi)
    e = gy(xi,yi,zi)
    f = gz(xi,yi,zi)
    g = hx(xi,yi,zi)
    h = hy(xi,yi,zi)
    i = hz(xi,yi,zi)
    Jcomplete = [a b c;d e f;g h i]
    return det(Jcomplete)
end

#test functions for R2
function f(x,y)
    fn = 2*x*(y&3)
    return fn
end

function g(x,y)
    fn = y*(x^2)
    return fn
end

function fx(x,y)
    fn = 2*(y^3)
    return fn
end

function fy(x,y)
    fn = 6*x*(y^2)
    return fn
end

function gx(x,y)
    fn = 2*x*y
    return fn
end

function gy(x,y)
    fn = x^2
    return fn
end

jacobian_eval_R2(fx,fy,gx,gy,1.,1.)

#test functions for R3
function f3(x,y,z)
    fn = 2*x*(y^3)*z
    return fn
end

function g3(x,y,z)
    fn = y*(x^2)*(z^4)
    return fn
end

function h3(x,y,z)
    fn = (y^4)*(z^3)
    return fn
end

function fx3(x,y,z)
    fn = 2*(y^3)*z
    return fn
end

function fy3(x,y,z)
    fn = 6*x*(y^2)*z
    return fn
end

function fz3(x,y,z)
    fn = 2*x*(y^3)
    return fn
end

function gx3(x,y,z)
    fn = 2*x*y*(z^4)
    return fn
end

function gy3(x,y,z)
    fn = x^2*(z^4)
    return fn
end

function gz3(x,y,z)
    fn = 4*(z^3)*y*(x^2)
    return fn
end

function hx3(x,y,z)
    fn = 0
    return fn
end

function hy3(x,y,z)
    fn = 4*(y^3)*(z^3)
    return fn
end

function hz3(x,y,z)
    fn = 3*(z^2)*(y^4)
    return fn
end

jacobian_eval_R3(fx3,fy3,fz3,hx3,gy3,gz3,hx3,hy3,hz3,1.,1.,1.)

#Wronskian Matrix evaluation at one point in R2.
#need to define f, f', g, g', and point (x,y).
function wronskian_eval_R2(f,fprime,g,gprime,x,y) #wronskian for a system of two equations.
    a = f(x,y)
    b = fprime(x,y)
    c = g(x,y)
    d = gprime(x,y)
    w = [a b;c d]
    if det(w) == 0
        print("linearly dependent")
    else
        print("linearly independent")
    end
    return det(w)
end


function f(x,y)
    fn = x^2
    return fn
end

function fprime(x,y)
    fn = 2*x*(y^3)
    return fn
end

function g(x,y)
    fn = x^3
    return fn
end

function gprime(x,y)
    fn = 3*(x^2)*y
    return fn
end

wronskian_eval_R2(f,fprime,g,gprime, 1., 1.) #this is the Wronskian only at one point.
