                                % n: number of neurons
global n;
n = 5;
                                % g: matrix of coupling conductances (n,n)
global g;
g = 2*rand(n)-1;
g(logical(eye(n))) = 0;
                                % tau: matrix of time delays (n,n)
tau = rand(n);
tau = .5*(tau + transpose(tau));

                                % lags: flattened list of time delays
lags = reshape(tau,[1,n^2]);

                                % vinf: equilibrium membrane potential
                                % winf: equilibrium recovery variable
                                % solve for equilibrium, select only real soln
global vinf;
vinf = roots([8 0 6 21]);
vinf = vinf(imag(vinf)==0);
global winf;
winf = vinf - vinf^3/3;
                                % y0: vector of equilibrium values
y0 = zeros(2*n,1);
y0(1:n) = vinf;
y0(n+1:2*n) = winf;

                                % sol: solution of dde
sol = dde23(@fndc,lags,y0,[0,60]);
plot(sol.x,sol.y(1:n,:));
xlabel('t (ms)')
ylabel('v (mV)')
saveas(gcf,'fndc.png')

                                % i_input: input current
                                % t: present time
function i = i_input(t)
  global n;
  i = zeros(n,1);
  if t>=1 && t<=2
    i(1) = 1;
  end
end

                           % fndc: (f)itzhugh-(n)agumo with (d)elayed (c)oupling
                           % t: present time
                           % y: present values of state variables (2n,1)
                           % z: past values of state variables (2n,n^2)
function dydt = fndc(t,y,z)
  global n;
  global g;
                                % v: present membrane potentials (n,1)
  v = y(1:n);
                                % w: present recovery variables (n,1)
  w = y(n+1:2*n);
                                % v_lag: past values of v (n,n)
                                % v_lag(i,j) = v_j(t - tau_{ij})
                                % k: 'flat' row-major index of (i,j)
  v_lag = zeros(n);
  for i = 1:n
    for j = 1:n
      k = n*(i-1) + j;
      v_lag(i,j) = z(j,k);
    end
  end
                                % i_network: current from network connections
  i_network = zeros(n,1);
  for i = 1:n
    i_network = i_network + g(:,i).*(v - v_lag(:,i));
  end

  dydt = [v - v.^3/3 - w - i_network + i_input(t)
          0.08*(v + 0.7 - 0.8*w)];
end
