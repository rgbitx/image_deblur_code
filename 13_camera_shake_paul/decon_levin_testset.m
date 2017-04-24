set_data_locations;

for imNum = 1:4
   for kerNum = 1:8
      imFile = ['im0' num2str(imNum) '_ker0' num2str(kerNum) '.mat'];
      test_blind_deconv_script;
   end
end