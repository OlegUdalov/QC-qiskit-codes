clear all
close all


theta_1 = 0.34;
phi_1 = 0.65;
lambda_1 = 1.4;

theta_2 = 2.34;
phi_2 = 3.1;
lambda_2 = 0.47;


%CZ_gate()*gates_product(U3_gate(theta_1, phi_1, lambda_1), U3_gate(theta_2, phi_2, lambda_2))


theta_3 = 0.65;
phi_3 = 0.83;
lambda_3 = 2.1;

theta_4 = 1.34;
phi_4 = 1.1;
lambda_4 = 0.17;

params_1 = [theta_1, phi_1, lambda_1, theta_2, phi_2, lambda_2, theta_3, phi_3, lambda_3, theta_4, phi_4, lambda_4];

init_circuit = single_cz_circ(params_1);

for i_try = 1:1
    options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',3000);
    %IterRec = @(params) (circ_diff(init_circuit, final_circuit(params)));
    IterRec = @(params) (circ_diff(init_circuit, double_cz_circ(params)));
    for i = 1:18
        params0(i)=pi*rand(1);
    end
    %circ_diff(init_circuit, final_circuit(params0))
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = -pi*ones(1,12);
    ub = pi*ones(1,12);
    nonlcon = [];
    params_fin=fmincon(IterRec,params0,A,b,Aeq,beq,lb,ub,nonlcon,options);

    %final_circuit(params_fin)
    %init_circuit
    %dd(i_try) = circ_diff(init_circuit, final_circuit(params_fin));
    dd(i_try) = (circ_diff(init_circuit, double_cz_circ(params_fin)));
%     for i = 1:6
%         fprintf([num2str(params_fin(i)),'\n'])
%     end
%     
%     fprintf([num2str(params_fin(7)), ' ', num2str(theta_1),'\n'])
%     fprintf([num2str(params_fin(8)), ' ', num2str(phi_1),'\n'])
%     fprintf([num2str(params_fin(9)), ' ', num2str(lambda_1),'\n'])
%     fprintf([num2str(params_fin(10)), ' ', num2str(theta_2),'\n'])
%     fprintf([num2str(params_fin(11)), ' ', num2str(phi_2),'\n'])
%     fprintf([num2str(params_fin(12)), ' ', num2str(lambda_2),'\n'])
    
    
    params_fin
    
    double_cz_circ(params_fin) 

end

dd

init_circuit

function diff = circ_diff_1(theta,phi,lambda, param)
    param_1(1) = theta;
    param_1(2) = phi;
    param_1(3) = lambda;
    param_1(4) = param(1);
    param_1(5) = param(2);
    param_1(6) = param(3);
    circ_1 = circ_sx(param_1);
    param_2(1) = param(4);
    param_2(2) = param(5);
    param_2(3) = param(6);
    param_2(4) = param(7);
    param_2(5) = param(8);
    param_2(6) = param(9);
    param_2(7) = param(10);
    param_2(8) = param(11);
    param_2(9) = param(12);
    circ_2 = circ_no_sx(param_2);
    diff = 0;
    for i = 1:4
        for j = 1:4
            diff = diff + (abs((circ_1(i,j)-circ_2(i,j))))^2;
        end
    end
    diff = sqrt(diff/16);
end

function circ = circ_sx(params)
    circ = gates_product(U3_gate(params(1), params(2), params(3)), U3_gate(params(4), params(5), params(6))) * CZ_gate() *  gates_product(I_gate(), sx_gate());
end

function circ = circ_no_sx(params)
    circ =  gates_product(U3_gate(params(1), params(2), params(3)), U3_gate(params(4), params(5), params(6))) * CZ_gate() *  gates_product(I_gate(), U3_gate(params(7), params(8), params(9)));
end

function circ = double_cz_circ(params)
    circ =  gates_product(U3_gate(params(1), params(2), params(3)), U3_gate(params(4), params(5), params(6))) * CZ_gate() *  gates_product(U3_gate(params(7), params(8), params(9)), U3_gate(params(10), params(11), params(12))) * CZ_gate() *  gates_product(U3_gate(params(13), params(14), params(15)), U3_gate(params(16), params(17), params(18)));
end

function circ = single_cz_circ(params)
    circ =  gates_product(U3_gate(params(1), params(2), params(3)), U3_gate(params(4), params(5), params(6))) * CZ_gate() *  gates_product(U3_gate(params(7), params(8), params(9)), U3_gate(params(10), params(11), params(12)));
end

function diff = circ_diff(circ_1, circ_2)
    diff = 0;
    for i = 1:4
        for j = 1:4
            diff = diff + (abs((circ_1(i,j)-circ_2(i,j))))^2;
        end
    end
    diff = sqrt(diff/16);
end


function circuit = final_circuit(params)
    circuit = gates_product(U3_gate(params(1), params(2), params(3)), U3_gate(params(4), params(5), params(6))) * CZ_gate() * gates_product(U3_gate(params(7), params(8), params(9)), U3_gate(params(10), params(11), params(12))) * CZ_gate() * gates_product(U3_gate(params(13), params(14), params(15)), U3_gate(params(16), params(17), params(18)));
end


function gate = I_gate()
    gate = zeros(2,2);
    gate(1,1) = 1;
    gate(2,2) = 1;
end

function DoubleQgate = CZ_gate()
    DoubleQgate = zeros(4,4);
    DoubleQgate(1,1) = 1;
    DoubleQgate(2,2) = 1;
    DoubleQgate(3,3) = 1;
    DoubleQgate(4,4) = -1;
end


function gate = U3_gate(theta_ang, phi_ang, lambda_ang)
    gate(1,1) = cos(theta_ang / 2);
    gate(1,2) = -sin(theta_ang / 2) * exp(1i * lambda_ang);
    gate(2,1) = sin(theta_ang / 2) * exp(1i * phi_ang);
    gate(2,2) = cos(theta_ang / 2) * exp(1i * (lambda_ang + phi_ang));
end

function gate = sx_gate()
    gate(1,1) = 1;
    gate(1,2) = 1i;
    gate(2,1) = 1i;
    gate(2,2) = 1;
end

function DoubleQgate = gates_product(G1, G2)
    DoubleQgate(1,1) = G1(1,1)*G2(1,1);
    DoubleQgate(1,2) = G1(1,1)*G2(1,2);
    DoubleQgate(2,1) = G1(1,1)*G2(2,1);
    DoubleQgate(2,2) = G1(1,1)*G2(2,2);
    
    DoubleQgate(1,3) = G1(1,2)*G2(1,1);
    DoubleQgate(1,4) = G1(1,2)*G2(1,2);
    DoubleQgate(2,3) = G1(1,2)*G2(2,1);
    DoubleQgate(2,4) = G1(1,2)*G2(2,2);
    
    
    DoubleQgate(3,1) = G1(2,1)*G2(1,1);
    DoubleQgate(3,2) = G1(2,1)*G2(1,2);
    DoubleQgate(4,1) = G1(2,1)*G2(2,1);
    DoubleQgate(4,2) = G1(2,1)*G2(2,2);
    
    
    DoubleQgate(3,3) = G1(2,2)*G2(1,1);
    DoubleQgate(3,4) = G1(2,2)*G2(1,2);
    DoubleQgate(4,3) = G1(2,2)*G2(2,1);
    DoubleQgate(4,4) = G1(2,2)*G2(2,2);
end