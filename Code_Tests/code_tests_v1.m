PLt = 2;
NDt = 4;
bt = 2 - PLt:1:0;
at = 2 - NDt:1:0;
%disp(a)
%disp(b)

av_datat = cell(272,1);
k = 2;
av_datat{k} = complex(600, 100);
av_datatk = av_datat{k};
%disp(av_datat)
%disp(length(av_datat))

%velocityt = zeros(length(av_datat)-1, length(at));
%disp(velocity)#

%nt = 1:length(av_datat)-2;
%disp(nt)
%jt = 1:length(bt);
%disp(jt)
%disp(bt(jt))

VD = zeros(length(av_datatk)-2,length(at));

%Matrix with length of av_data cell rows and 2 columns.
v = zeros(length(av_datatk)-2,2);

%Matrix with length of av_data cell rows and mm (no. of gates) columns.
V_dop = zeros(length(av_datatk)-2,5);

for n = 1:length(av_datatk)-2
    m = 1;

    %k = 1, then k = 2, k = 3. Loops three times.
    for k = 1:length(at)
        %Index n,k in matrix VD and set it to 0.
        VD(n,k) = 0;

        %Only loops once, length of b is 1.
        for j = 1:length(bt)

            %This if statement is for the first iteration.
               if n==1 && bt(j) == 0

                   %Index j,k in matrix v and set it to output of the
                   %operation below.
                   disp(m-at(k))
                   v(j,k) = av_datatk(n+2-bt(j),...
                       m-at(k)).*conj(av_datatk(n+1-bt(j),m-at(k)));
                   VD(n,k) = VD(n,k)+ v(j,k);
                   disp("bt(j) is 0 lel")
                   disp(VD)
               elseif n==1 && bt(j) ==1
                   v(j,k) = av_datat(n-1+bt(j),...
                       m-at(k)).*conj(av_datat(n-1 + bt(j)+1,m-at(k)));
                   VD(n,k) = VD(n,k) + v(j,k);
                   disp("bt(j) is 1 lel")
                   disp(VD)
               else
                   v(j,k) = av_datat(n-bt(j)+1,...
                       m-at(k)).*conj(av_datat(n-bt(j),m-at(k)));
                   VD(n,k) = VD(n,k)+ v(j,k);
                   disp("bt(j) is neither 0 or 1 lel")
                   disp(VD)
               end
        end
        V_dop(n,m) = V_dop(n,m) + VD(n,k);
        disp(V_dop)
    end
end
