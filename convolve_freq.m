[input, fs] = wavread('voice.wav');
ir = wavread('ir.wav');

inputlen = length(input);
irlen = length(ir);

inputbuffered = zeros(irlen, 1);
inputbuffered(1:inputlen) = input;

% no windowing yet...
inputspec = fft(inputbuffered);
irspec = fft(ir);

outputspec = inputspec .* irspec;

output = ifft(outputspec);

divisor = max(abs(output));
divisor = 1.0/divisor; % use recip to avoid divides below

% normalize
output = output*divisor;

wavwrite(output,fs,'output.wav');