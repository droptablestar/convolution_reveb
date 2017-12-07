[input, fs] = wavread('continuous.wav');
ir = wavread('ir.wav');

inputlen = length(input);
irlen = length(ir);
outputlen = inputlen+irlen-1; % from Moore intro part 2, p. 53

output = zeros(outputlen, 1);

for i=1:inputlen,
	output(i:i+irlen-1) = output(i:i+irlen-1) .+ input(i)*ir;
end

divisor = max(abs(output));
divisor = 1.0/divisor; % use recip to avoid divides below

% normalize
output = output*divisor;

wavwrite(output,fs,'output.wav');