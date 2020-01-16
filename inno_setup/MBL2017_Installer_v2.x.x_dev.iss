; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{034BDD77-7FEF-4D05-9038-04FE4C4B168D}
AppName=MBL2017
AppVersion=2.0.0
AppPublisher=MBL2017 Software by Palaeoware
AppPublisherURL=https://github.com/palaeoware/MBL2017/
AppSupportURL=https://github.com/palaeoware/MBL2017/
AppUpdatesURL=https://github.com/palaeoware/MBL2017/
DefaultDirName={pf}\Palaeoware\MBL2017\v2.0.0
DefaultGroupName=Palaeoware\MBL2017
ArchitecturesAllowed=x64
ArchitecturesInstallIn64BitMode=x64
LicenseFile=.\bin\LICENSE.md
OutputDir=.\build
OutputBaseFilename=MBL2017Installer_v2.0.0_win_x64
Compression=lzma
SolidCompression=yes
ChangesAssociations=yes
AppCopyright=Copyright (C) 2016-2019 by Mark D. Sutton, Julia Sigwart, Russell J. Garwood, and Alan R.T. Spencer
AppContact=palaeoware@gmail.com
AppComments=MBL2017 (Phylogenetic tree and matrix simulation based on birth–death systems). More details are available in the readme for the program.
BackColor=$3e3e3e
BackColor2=$3e3e3e

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Files]
Source: ".\bin\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
; NOTE: Don't use "Flags: ignoreversion" on any shared system files

[Icons]
Name: "{group}\MBL2017 (v2.0.0)"; Filename: "{app}\MBL2017.exe";

[Run]
Filename: "{app}\MBL2017.exe"; Description: "{cm:LaunchProgram,MBL2017}"; Flags: nowait postinstall skipifsilent