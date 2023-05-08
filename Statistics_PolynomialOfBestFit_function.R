# Michael Gunnarson
# Summer 2022
# Written in R 4.2.0

poly_bf <- function(x,y) {
	# input X and Y column data.  algorithmically solve for statistically significant
	# polynomial of best fit (alpha = .05)  output x and y data for line of best fit, 
	# as well as coefficient matrix of the form y = a + bx + cx^2 + ...
	
	str0 = 'y~x+0'
	str1 = 'y~x'
	master_string = " + I(x^n)"



	# initialize
	i = 0

	p = 0
	alpha = 0.05

	repl = "n"

	THE_ANSWER <- list("placeholder1","placeholder2") 

	# a rotating window of strings, two variables wide
	# at the end of the loop, we will have gone one too far.  The fit information cannot, as far as I know, be
	# stored in sequential variables or an array of arrays.  If we keep a rolling window on the last two vectors, 
	# we can see the variable that pushed us too far and the one that is the line of best fit.
	# thus we export the smaller of the two strings and use it to re-solve for the line of best fit


	Caterpiller <- function(list,replace) {
	# takes 2 element long list and slides it to the left
	# crude circshift and replace

		list[1] = list[2]
		list[2] = replace
		return(list)
	}


	# general while loop for linear model of best fit up to polynomial of order 999


	while (p < alpha) {
		i = i+1 # i will be one too big when loop breaks

		if (i ==  1) {
			fit <- lm(str0) # works with string variables!!!
			p = summary(fit)$coefficients[i,4] # check value, Pr >t

			THE_ANSWER = Caterpiller(THE_ANSWER,str0)
			}

		else if (i == 9) {
			fit <- lm(str1)
			p = summary(fit)$coefficients[i,4] # check value, Pr >t
	
			THE_ANSWER = Caterpiller(THE_ANSWER,str1)
	
			# deal with string case n >10
			master_string = " + I(x^nn)"
			repl = "nn"
	
			# replace n in string with n variable
			str_conc = sub(repl, i, master_string) 
			
			str1 = paste(str1, str_conc) # continually concatenate strings each loop
			}
	
		else if (i == 99) {
			fit <- lm(str1)
			p = summary(fit)$coefficients[i,4] # check value, Pr >t
	
			THE_ANSWER = Caterpiller(THE_ANSWER,str1)
			
			# deal with string case n >100
			master_string = " + I(x^nnn)"
			repl = "nnn"
	
			# replace n in string with n variable
			str_conc = sub(repl, i, master_string) 
			
			str1 = paste(str1, str_conc) # continually concatenate strings each loop
			}
	
		else if (i == 1000) {
			break
			}
	
		else {
			fit <- lm(str1)
			p = summary(fit)$coefficients[i,4] # check value, Pr >t
			
			THE_ANSWER = Caterpiller(THE_ANSWER,str1)
	
			# replace n in string with n variable
			str_conc = sub(repl, i, master_string) 
			
			str1 = paste(str1, str_conc) # continually concatenate strings each loop
			}
	# 	print(i)
	#	print(summary(fit))
	#	print(p)
	#	print(THE_ANSWER[[2]])  # by inspection, this data set is indeed a polynomial of degree x^5
						# our loop works!
	}
	
	
	# I don't need to keep track of every single fit, especially since i can't keep fit data in a matrix or array.  On the other hand,
	# I can in fact store the strings in a list
	
	
	optimal_fit <- lm(THE_ANSWER[[1]])
	unoptimal_fit <- lm(THE_ANSWER[[2]])
	
	n = i-1 # if i = 7, polynomial of degree 6.  if i = 7, desired i is 6, degree is 5
	# n is number of terms
	
	deg = n-1
	
	
	coef = optimal_fit$coefficients
	yfit = 0
	xfit = seq(from = min(x), to = max(x), length.out = 50)
	
	for (j in 1:n){
		yfit = yfit + coef[j]*xfit^(j-1)
	}
		
	return(list(xfit,yfit,coef))
}
