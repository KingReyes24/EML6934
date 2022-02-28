function state = myode(t,X)

% state = 5*[
state = [
    -X(1)-X(2); % xDot       = -x - lambda
    -X(1)+X(2); % lambdaDot  = -x + lambda
    ];

end