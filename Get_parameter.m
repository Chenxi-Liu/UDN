function [m_l,counter] = Get_parameter(n)
a=1:n;
%
counter = 1;
%%
if n==1
    for m_1 = 0:n
        if a*[m_1]'==n
            m_l(:,:,counter) = [m_1];
            counter = counter + 1;
        end
    end
end
%%
if n==2
    for m_1 = 0 : n
        m_1
        for m_2 = 0 : n
            if a*[m_1,m_2]'==n
            m_l(:,:,counter) = [m_1,m_2];
            counter = counter + 1;
            end
        end
    end
end
%%
if n==3
    for m_1 = 0 : n
        m_1
        for m_2 = 0 : n
            for m_3 = 0 : n
                if a*[m_1,m_2,m_3]'==n
                    m_l(:,:,counter) = [m_1,m_2,m_3];
                    counter = counter + 1;
                end
            end
        end
    end
end
%%
if n==4
    for m_1 = 0 : n
        m_1
        for m_2 = 0 : n
            for m_3 = 0 : n
                for m_4 = 0 : n
                    if a*[m_1,m_2,m_3,m_4]'==n
                        m_l(:,:,counter) = [m_1,m_2,m_3,m_4];
                        counter = counter + 1;
                    end
                end
            end
        end
    end
end
%%
if n==5
    for m_1 = 0 : n
        m_1
        for m_2 = 0 : n
            for m_3 = 0 : n
                for m_4 = 0 : n
                    for m_5 = 0 : n
                        if a*[m_1,m_2,m_3,m_4,m_5]'==n
                            m_l(:,:,counter) = [m_1,m_2,m_3,m_4,m_5];
                            counter = counter + 1;
                        end
                    end
                end
            end
        end
    end
end
if n==6
    for m_1 = 0:n
        m_1
        for m_2 = 0:n
            for m_3 = 0:n
                for m_4 = 0:n
                    for m_5 = 0:n
                        for m_6 = 0:n
                            if a*[m_1,m_2,m_3,m_4,m_5,m_6]'==n
                                m_l(:,:,counter)=[m_1,m_2,m_3,m_4,m_5,m_6];
                                counter = counter + 1;
                            end                        
                        end
                    end
                end
            end
        end
    end
end
if n==7
    for m_1 = 0:n
        m_1
        for m_2 = 0:n
            for m_3 = 0:n
                for m_4 = 0:n
                    for m_5 = 0:n
                        for m_6 = 0:n
                            for m_7 = 0:n
                                if a*[m_1,m_2,m_3,m_4,m_5,m_6,m_7]'==n
                                    m_l(:,:,counter)=[m_1,m_2,m_3,m_4,m_5,m_6,m_7];
                                    counter = counter + 1;
                                end     
                            end
                        end
                    end
                end
            end
        end
    end
end
if n==8
   for m_1 = 0:n
       m_1
        for m_2 = 0:n
            for m_3 = 0:n
                for m_4 = 0:n
                    for m_5 = 0:n
                        for m_6 = 0:n
                            for m_7 = 0:n
                                for m_8 = 0:n
                                    if a*[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8]'==n
                                        m_l(:,:,counter)=[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8];
                                        counter = counter + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if n==9
   for m_1 = 0:n
       m_1
        for m_2 = 0:n
            for m_3 = 0:n
                for m_4 = 0:n
                    for m_5 = 0:n
                        for m_6 = 0:n
                            for m_7 = 0:n
                                for m_8 = 0:n
                                    for m_9 = 0 : n
                                        if a*[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9]'==n
                                            m_l(:,:,counter)=[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9];
                                            counter = counter + 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if n==10
   for m_1 = 0:n
       m_1
        for m_2 = 0:n
            for m_3 = 0:n
                for m_4 = 0:n
                    for m_5 = 0:n
                        for m_6 = 0:n
                            for m_7 = 0:n
                                for m_8 = 0:n
                                    for m_9 = 0 : n
                                        for m_10 = 0:n
                                            if a*[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10]'==n
                                                m_l(:,:,counter)=[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10];
                                                counter = counter + 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if n==11
   for m_1 = 0:n
       m_1
        for m_2 = 0:n
            for m_3 = 0:n
                for m_4 = 0:n
                    for m_5 = 0:n
                        for m_6 = 0:n
                            for m_7 = 0:n
                                for m_8 = 0:n
                                    for m_9 = 0 : n
                                        for m_10 = 0:n
                                            for m_11=0:n
                                                if a*[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10,m_11]'==n
                                                    m_l(:,:,counter)=[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10,m_11];
                                                    counter = counter + 1;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if n==12
   for m_1 = 0:n
       m_1
        for m_2 = 0:n
            for m_3 = 0:n
                for m_4 = 0:n
                    for m_5 = 0:n
                        for m_6 = 0:n
                            for m_7 = 0:n
                                for m_8 = 0:n
                                    for m_9 = 0 : n
                                        for m_10 = 0:n
                                            for m_11=0:n
                                                for m_12=0:n
                                                    if a*[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10,m_11,m_12]'==n
                                                        m_l(:,:,counter)=[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10,m_11,m_12];
                                                        counter = counter + 1;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if n==13
   for m_1 = 0:n
       m_1
        for m_2 = 0:n
            for m_3 = 0:n
                for m_4 = 0:n
                    for m_5 = 0:n
                        for m_6 = 0:n
                            for m_7 = 0:n
                                for m_8 = 0:n
                                    for m_9 = 0 : n
                                        for m_10 = 0:n
                                            for m_11=0:n
                                                for m_12=0:n
                                                    for m_13 = 0:n
                                                        if a*[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10,m_11,m_12,m_13]'==n
                                                            m_l(:,:,counter)=[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10,m_11,m_12,m_13];
                                                            counter = counter + 1;
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if n==14
   for m_1 = 0:n
       m_1
        for m_2 = 0:n
            for m_3 = 0:n
                for m_4 = 0:n
                    for m_5 = 0:n
                        for m_6 = 0:n
                            for m_7 = 0:n
                                for m_8 = 0:n
                                    for m_9 = 0 : n
                                        for m_10 = 0:n
                                            for m_11=0:n
                                                for m_12=0:n
                                                    for m_13 = 0:n
                                                        for m_14=0:n
                                                            if a*[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10,m_11,m_12,m_13,m_14]'==n
                                                                m_l(:,:,counter)=[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10,m_11,m_12,m_13,m_14];
                                                                counter = counter + 1;
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if n==15
   for m_1 = 0:n
       m_1
        for m_2 = 0:n
            for m_3 = 0:n
                for m_4 = 0:n
                    for m_5 = 0:n
                        for m_6 = 0:n
                            for m_7 = 0:n
                                for m_8 = 0:n
                                    for m_9 = 0 : n
                                        for m_10 = 0:n
                                            for m_11=0:n
                                                for m_12=0:n
                                                    for m_13 = 0:n
                                                        for m_14=0:n
                                                            for m_15 = 0:n
                                                                if a*[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10,m_11,m_12,m_13,m_14,m_15]'==n
                                                                    m_l(:,:,counter)=[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10,m_11,m_12,m_13,m_14,m_15];
                                                                    counter = counter + 1;
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if n==16
   for m_1 = 0:n
       m_1
        for m_2 = 0:n
            for m_3 = 0:n
                for m_4 = 0:n
                    for m_5 = 0:n
                        for m_6 = 0:n
                            for m_7 = 0:n
                                for m_8 = 0:n
                                    for m_9 = 0 : n
                                        for m_10 = 0:n
                                            for m_11=0:n
                                                for m_12=0:n
                                                    for m_13 = 0:n
                                                        for m_14=0:n
                                                            for m_15 = 0:n
                                                                for m_16=0:n
                                                                    if a*[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10,m_11,m_12,m_13,m_14,m_15,m_16]'==n
                                                                        m_l(:,:,counter)=[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10,m_11,m_12,m_13,m_14,m_15,m_16];
                                                                        counter = counter + 1;
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
counter = counter-1;