%����Mackay�Ĺ��취1A����У�����
function [H]=getH1(m,n)
row_flag(1:m)=0;
h=zeros(m,n);
bits_per_col=4;                                                            %ʹÿ���������4��1������Ϊ4
for i=1:n
    a=randperm(m);
    for j=1:bits_per_col
        h(a(j),i)=1;
      row_flag(a(j))= row_flag(a(j))+1;
    end
end                                                                        %ÿ�ж���2��1��ÿ��1��λ��������Ĳ�ͬ��
max_ones_per_row=ceil(n*bits_per_col/m);                                   %��ÿ��1�������������ص����ֵ
for i=1:m                                                                  
    if row_flag(i)==0
        for k=1:3   %����
            j=unidrnd(m);                                                  %j=1~m�������ֵ
            while  h(i,j)==1
                j=unidrnd(m);
            end
            h(i,j)=1;
            row_flag(i)=row_flag(i)+1;
        end
    end
    if  row_flag(i)==1
         j=unidrnd(m);
            while  h(i,j)==1
                j=unidrnd(m);
            end
            h(i,j)=1;
            row_flag(i)=row_flag(i)+1;
    end
end
          
for i=1:m                                                                  %���������Ϸ�ɢ1��λ��,ʹ���طֲ���������
    j=1;
    a=randperm(n);
    while row_flag(i)>max_ones_per_row                                     %����������ش������ص����ֵ,����д���
        if h(i,a(j))==1                                                    %���ѡ��ĳһ������Ϊ1����������,�������ϵ�1��ɢ����������
            newrow=unidrnd(m);                                             %������Ҹ������ʺϷ���1(����С�����ֵ�Ҹ�λ��Ϊ0)����
            k=0;
            while (row_flag(newrow)>=max_ones_per_row|h(newrow,a(j))==1)&k<m
                newrow=unidrnd(m);
                k=k+1;
            end
            if h(newrow,a(j))==0                                           %���������е�1�ŵ��ҵ�������,�������е��б�־����Ӧ�Ĵ���
                h(newrow,a(j))=1;
                row_flag(newrow)=row_flag(newrow)+1;
                h(i,a(j))=0;
                row_flag(i)=row_flag(i)-1;
            end
        end
        j=j+1;
    end
end

for loop=1:100                                     %����ɾ���̻�4
    success=1;
    for r=1:m
        ones_position=find(h(r,:)==1);             %����r��Ϊ1��Ԫ�ص�λ���ҵ�����Ϊones_position
        ones_count=length(ones_position);
        for i=[1:r-1,r+1:m]
            common=0;
            for j=1:ones_count
                if h(i,ones_position(j))==1
                    common= common+1;
                    if  common==1
                        a=ones_position(j);         %���еĵ�һ����ͬ1Ԫ��
                    end
                end
                if common==2
                    success=0;
                    common=common-1;
                    if (round(rand)==0)            %�漴��������ǰ����л��Ǻ������
                        b=a;                       %�����������,����ǰ�����ϵĵڶ�����ͬ��1Ԫ��
                        a=ones_position(j);
                    else b=ones_position(j);       %����ǰ�����,����������ϵĵڶ�����ͬ��1Ԫ��
                    end
                    h(i,b)=3;                      %����1��Ϊ3��ʹ���Ժ�ĳ���ɾ���в��ø�ֵ
                    newrow=unidrnd(m);
                    iteration=0;
                    while h(newrow,b)~=0 & iteration<5 %����5���ڴ������������������0
                        newrow=unidrnd(m);
                        iteration=iteration+1;
                    end
                    if iteration>=5                   %����5����������ҷ�1��0��3
                        while h(newrow,b)==1
                           newrow=unidrnd(m);
                        end
                    end
                    h(newrow,b)=1;                    %���������ҵ���0��3��Ϊ1
                end
            end
        end
    end
    if success                                        %�������ѭ���Ѳ����ڶ̻�4(seccess=1),�����ѭ��loop
        break
    end
end
h=h==1;                                               %��0���ʣ���3,���õ�����õ�У�����H
H=h;
%save H;