#include "mainwindow.h"
#include "version.h"
#include "darkstyletheme.h"

#include <QApplication>
#include <QSplashScreen>
#include <QString>
#include <QStyle>
#include <QDesktopWidget>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

   //Close on last window close
   a.setQuitOnLastWindowClosed(true);

   //Style program with our dark style
   a.setStyle(new DarkStyleTheme);

   QPixmap splashPixmap(":/palaeoware_logo_square.png");
   QSplashScreen splash(splashPixmap,Qt::WindowStaysOnTopHint);
   splash.show();
   QString version;
   version.sprintf("%d.%d.%d",MAJORVERSION,MINORVERSION,PATCHVERSION);
   splash.showMessage("<font><b>" + QString(PRODUCTNAME) + " - v" + version + "</b></font>",Qt::AlignHCenter,Qt::white);
   a.processEvents();


    MainWindow w;
    w.show();

    return a.exec();
}
