function result = AlMin_MC()

Sparsity_range = 0.1:0.025:0.975;
N_range = 20:10:200;
%epsilon = epsilon_list(epsilon_index);
%sigma = (I^2*sqrt(256*T*log(2/delta)))/epsilon;
result = zeros(length(Sparsity_range),length(N_range));
for m = 1:length(Sparsity_range)
    for n = 1:length(N_range)
        M = 200;
        L = 10;
        N = N_range(n);
        A = randn(M,L);
        X = randn(L,N);
        D = A*X;
        [MM,NN] = size(D); %返回data数据文件里的矩阵大小
        rho = Sparsity_range(m);
        Omega = rand(MM,NN)<=rho; % support of observation
        D_omega = Omega.*D; % measurements are made
        D_omega1= D_omega';
        Omega1 =Omega';
        PP=rand(1);
        P=PP(1,1);
        RR=rank(D);
        T=10;
        I = maxl2norm(D,Omega);
        delta=2.2251e-308;
        %delta=10^(-6);
        temp=sqrt(2*log(2/(delta)))/(2*log(1/(delta)));
        sigma = 2*I*T*sqrt(2*log(2/(delta*T)))/(2*log(1/(delta)));
        epsilon_list = [1];
        x1 = zeros(M,RR);
        x2 = zeros(N,RR);
        v = zeros(RR,NN);
        u = zeros(RR,MM);
        for t=1:T  %循环次数
            if t==1
                %U = randi([0,5],MM,RR);
                temp=D_omega.*P;
                [U,S1,V1]=svds(temp,RR);
            end

            for j=1:NN %分别算V的每一列
                y=D_omega(:,j);        
                for i=1:RR
                    x1(:,i)=Omega(:,j).*U(:,i);%将采样矩阵K的第j列跟U0对应行的每个元素相乘
                end
                v(:,j)=pinv(x1)*y; %算v的第j列
            end
            V=v';
            for j=1:MM %分别算U的每一行
                y=D_omega1(:,j);
                for i=1:RR
                    x2(:,i)=Omega1(:,j).*V(:,i);
                end
                u(:,j)=pinv(x2)*y; %算A'的第j列
            end
            U=u';V=V';
            U=U+normrnd(0,sigma^2*T);
            M0=U*V;  
            result(m,n) = norm((D-M0),'fro')/(M*N);
            disp(['Sparsity_range:' num2str(Sparsity_range(m)) 'N_range:'...
                num2str(N_range(n)) 'result:' num2str(result(m,n))]);
            %p2(t)=rmse(D,M0);        
        end
    end
end
save('C:\Users\DELL\Desktop\Paper use code\_result\MC_AlMin_Phase_Transiton','result');
end
%save('C:\Users\DELL\Desktop\Paper use code\_result\MC_AlMin_process','result');
%result(epsilon_index) = p2(T);
%plot(1:T,result);
% 20*log10(p(10))
%q=[1:T];
% plot(epsilon_list, result, '-bs','MarkerEdgeColor','r','MarkerFaceColor','r', 'MarkerSize', 10);
% xlabel('Epsilon per user');
% ylabel('RMSE');
% grid on;
%pic=semilogy(q,p2);