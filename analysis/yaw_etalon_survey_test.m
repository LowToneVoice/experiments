// target_obj:yaw_etalon
// command:setPosition
system startLogging;

yaw_etalon setPosition:6.480000;
niki startDetection;
gate waitWhileKP:1500;
niki stopDetection;

yaw_etalon setPosition:6.465000;
niki startDetection;
gate waitWhileKP:1500;
niki stopDetection;

yaw_etalon setPosition:6.480000;
niki startDetection;
gate waitWhileKP:1500;
niki stopDetection;

yaw_etalon setPosition:6.475000;
niki startDetection;
gate waitWhileKP:1500;
niki stopDetection;

yaw_etalon setPosition:6.480000;
niki startDetection;
gate waitWhileKP:1500;
niki stopDetection;

yaw_etalon setPosition:6.485000;
niki startDetection;
gate waitWhileKP:1500;
niki stopDetection;

system stopLogging;
gate GatenetLogRename:;
