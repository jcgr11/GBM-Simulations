# Monte Carlo simulation for stock price dynamics using a recursive form of the Geometric Brownian Motion SDE
* The purpose of this project is to quantify and vizualize the uncertainty and randomness inherent in stock price movement.

## References:
* Geometric brownian motion. (n.d.). https://www.quantstart.com/articles/Geometric-Brownian-Motion/
* Lewinson, E. (2022). Python for finance cookbook: Over 80 powerful recipes for effective financial data analysis. Packt Publishing.
* K. (2013). Simulating brownian motion (bm) and Geometric Brownian motion (GBM). http://www.columbia.edu/~ks20/4404-Sigman/4404-Notes-sim-BM.pdf
* Rouah, F. (n.d.). Euler and Milstein discretization - frouah.com. http://www.frouah.com/finance%20notes/Euler%20and%20Milstein%20Discretization.pdf
* MIT OpenCourseWare. (n.d.). Lecture 21: Stochastic differential equations: Topics in mathematics with applications in finance: Mathematics. MIT OpenCourseWare. Retrieved April 4, 2023, from https://ocw.mit.edu/courses/18-s096-topics-in-mathematics-with-applications-in-finance-fall-2013/resources/lecture-21-stochastic-differential-equations/ 

## A breif description of Geometric brownian motion, the derived recursive form used in this model for estimating geometric borwnian motion in stock price dynamics:
### Geometric brownian motion:

Geometric Brownian Motion is a continuous time stochastic processs used to describe the stochastic movement of stock prices. 

Stock prices are considered to be a stochastic process beause they are subject to random fluctuations that are influenced by large number of uncontrollable factors (i.e. economic new, company performance, and investor sentiment). 

As with any quantitive model used for informing investment decisions there are limitations. The Geometric Brownian Motion Model assumes that the volatility of the price process is constant over time, and that returns are normally distributed across the observed time interval. 

Even after considering the limitations of the model, it continues to provide a well defined picture of what to expect in terms of how a particular stock may behave in the financial markets.

### The GBM model is expressed by the following stochastic differential equation:
<p align="center">
$$dS(t) = μS(t)dt + σS(t)dW(t)$$
</p>

where:
* $S(t)$ is the price of the asset at time t
* $μ$, aka the drift coefficient, is the expected rate of return on the asset
* $σ$, aka the diffusion, is the volatility of the asset
* $W(t)$ is a Wiener process, which is a mathematical representation of Brownian motion

In this case W(t) is simulated using using a simple Monte Carlo simulation by generating a series normally distrubted random numbers which are used to update the state of the stochastic process at each time step, although not an exact weiner process the simulation will have similar properties such as a continuous and unpredicatable path at each time step.

### The derivation of the recursive form of the GBM SDE (shown below) using Euler-Maruyama's Discretization:
<p align="center">
$$S(t+dt)=S(t)\cdot\exp\left((\mu-\frac{1}{2}\sigma^2)dt+\sigma\sqrt{dt}\cdot Z\right)$$
</p>

where:
* $S(t)$ is the stock price at time $t$
* $dt$ is the time step size
* $μ$ is the drift or expected stock return
* $σ$ is the volatility or diffusion
* $Z$ is a standard normal random variable
* $exp$ is the exponential function

and:
* $(\mu-\frac{1}{2}\sigma^2)dt$ represents the expected change in the log price over the time step

and:
* $\sigma\sqrt{dt}\cdot Z$  represents the random component of the log price change.

This recursive form of the GBM SDE is derived by discretizing the SDE using the Euler-Maruyama method, which approximates the stochastic differential equation by a difference equation.
* Note: A difference equation is a mathematical equation that describes a sequence of values in terms of the values that precede it.

#### Euler discretization to discretize the GBM SDE: $dS(t) = μS(t)dt + σS(t)dW(t)$

* Start by defining the Wiener process increment $dW(t)$ as a standard normal random variable times the square root of the time step $dt$:
$$dW(t) = \sqrt{dt}\cdot Z$$
where Z is a standard normal random variable.

* Then, approximate the change in $S(t)$ over the time step dt as:
$$dS(t) = μS(t)dt + σS(t)\sqrt{dt}\cdot Z$$$
This is the same as the original SDE, but with $dW(t)$ replaced by our approximation from step 1.

Rearrange the terms in the above equation to isolate $S(t + dt)$:

* Start with the Original SDE where $\sqrt{dt} Z$ is used to estimate dW(t)$: $$dS(t) = \mu S(t) dt + \sigma S(t) \sqrt{dt} Z$$
* Divide both side by $S(t)$: $$\frac{dS(t)}{S(t)} = \mu dt + \sigma \sqrt{dt} Z$$
* Subtract both sides by $\mu dt$: $$\frac{dS(t)}{S(t)} - \mu dt = \sigma \sqrt{dt} Z$$
* $d\ln(S(t)) - \mu dt = \sigma \sqrt{dt} Z$
* $\ln(S(t+dt)) - \ln(S(t)) = \mu dt + \sigma \sqrt{dt} Z$
* $\ln(S(t+dt)) = \ln(S(t)) + \mu dt + \sigma \sqrt{dt} Z$
* $S(t+dt) = S(t) + \mu S(t) dt + \sigma S(t) \sqrt{dt} Z$

where $ln()$ is the natural logarithm function.

* Now, notice that the term $\mu S(t) dt$ can be written as:

$$\mu S(t) dt = \mu S(t) dt - \frac{1}{2}\sigma^2 S(t) dt + \frac{1}{2}\sigma^2 S(t) dt$$

$$= (\mu - \frac{1}{2}\sigma^2) S(t) dt + \frac{1}{2}\sigma^2 S(t) dt$$

* So, the discretized equation can be rewritten as:

$$S(t+dt) = S(t) + (\mu - \frac{1}{2}\sigma^2) S(t) dt + \sigma S(t) \sqrt{dt} Z$$

* Finally, this expression can be simplified by factoring out $S(t)$:

$$S(t+dt) = S(t) \left(1 + (\mu - \frac{1}{2}\sigma^2) dt + \sigma \sqrt{dt} Z \right)$$
* which is equivalent to:

$$S(t+dt) = S(t) \exp \left( (\mu - \frac{1}{2}\sigma^2) dt + \sigma \sqrt{dt} Z \right)$$

The final equation is the discrete approximation to the original GBM SDE, using the Euler-Maruyama method. It tells us that the predicted stock price after a small time interval $dt$, $S(t + dt)$, can be obtained by multiplying the current stock price $S(t)$ by a random factor $\exp \left( (\mu - \frac{1}{2}\sigma^2) dt + \sigma \sqrt{dt} Z \right)$, where $Z$ is a standard normal random variable. This random factor captures the stochastic fluctuations in the stock price due to the Wiener process $dW(t)$.

# The code:

## Setup:
* Import necessary packages (see requirements file).
* Define the security for which you would like to simulate price dynamics for. 
* Choose the time frame from which you will train and test the Monte Carlo Simulation model.
* Download you data using the yf.download() function.
* specify the Adj Close column to account for an income accrual events in the price of the security (dividends, short term and long term capital gains, corp actions, etc.). adj_close = x['Adj Close']
* Calculate the daily percentage change in adjusted close prices and drop the first row which contains null values. returns = adj_close.pct_change().dropna()

## Setting up the model:
* Split your returns dataframe into your training set and test set
* Set up your Monte Carlo simulation parametes: Including the length of the test period (T), the number of time steps (N), the initial stock price (S_0), the number of simulations (N_SIM), and the mean (mu) and standard deviation (sigma) of the returns during the training period.
### The simulate_gbm function:
* This function simulates the Geometric Brownian Motion (GBM) of the selcted security using a Monte Carlo simulation. The function input parameters are :
* s_0: the initial stock price
* mu: the expected return of the stock
* sigma: the volatility of the stock
* n_sims: the number of simulations to run
* T: the time horizon
* N: the number of time steps to divide the horizon into

The function calculates the time step dt by dividing T by N. It then generates an array of random normal variables dW with a scale of the square root of dt using the NumPy np.random.normal function.

Next, the simulate_gbm function uses NumPy's np.cumsum function to compute the cumulative sum of dW along the second axis to obtain an array W of Brownian motion values. It also creates a matrix time_steps consisting of repeated copies of the time steps.

Finally, the function computes the simulated stock prices S_t at each time step by applying the GBM formula using the np.exp and np.insert functions to insert the initial stock price s_0 at the beginning of each row of the array. The function returns the simulated stock prices S_t.

### Preparing and Vizualizing your results against to original mean return:
* LAST_TRAIN_DATE = train.index[-1].date() - sets the variable LAST_TRAIN_DATE to the date of the last entry in the train DataFrame's index. train is a DataFrame containing historical data up until the last date of the training set.
* FIRST_TEST_DATE = test.index[0].date() - sets the variable FIRST_TEST_DATE to the date of the first entry in the test DataFrame's index. test is a DataFrame containing historical data starting from the first date of the testing set.
* LAST_TEST_DATE = test.index[-1].date() - sets the variable LAST_TEST_DATE to the date of the last entry in the test DataFrame's index.
* PLOT_TITLE = (f'{ticker} Simulation 'f'({FIRST_TEST_DATE}:{LAST_TEST_DATE})') - creates a string that will be used as the title of the plot. 
* selected_indices selects the date range for which we have both actual stock price data and simulated data.
* gbm_simulations_df creates a pandas DataFrame of the GBM simulation results, with the date range from selected_indices as the index.
* The plot function is then used to plot the DataFrame gbm_simulations_df with a specified figsize and an alpha of 0.2 (to make the simulation lines more transparent), and with no legend (legend=False).
* line_1 and line_2 are created using the plot function, with line_1 representing the simulated mean of the GBM simulations, and line_2 representing the actual stock price path.
* The set_title function sets the plot title to PLOT_TITLE.
* The legend function adds a legend to the plot, with line_1 labeled as 'Simulated Mean' and line_2 labeled as 'Actual Price Path'.
* The plot is then displayed using plt.show().
