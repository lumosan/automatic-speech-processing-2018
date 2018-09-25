% FUNCTION_MAPPER Initializes function handles so that different implementations of a function can be easily interchanged
% Comment/uncomment one line of every pair to choose between your implementation and the compiled solution


    hw1_ammod = str2func('sol_ammod');
%     hw1_ammod = str2func('my_ammod');  % uncomment to use your solution

    hw1_amdemod = str2func('sol_amdemod');
%     hw1_amdemod = str2func('my_amdemod'); % uncomment to use your solution