function ableTemp(~,~,f)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
fh = guidata(f);
if(fh.needTemp.Value)
    fh.tempString.Visible = 'on';
    fh.tempSetting.Visible = 'on';
    fh.needRefine.Visible = 'on';
    fh.needGrow.Visible = 'on';
else
    fh.tempString.Visible = 'off';
    fh.tempSetting.Visible = 'off';
    fh.needRefine.Visible = 'off';
    fh.needGrow.Visible = 'off';
end
end

