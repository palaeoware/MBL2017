#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "simulation.h"
#include "lineage.h"
#include "qcustomplot.h"
#include <QHash>
#include <QFileDialog>
#include <QPen>
#include <QMessageBox>


/////////////////////////////////////////////////////
//MainWindow object - handle GUI and user actions
/////////////////////////////////////////////////////


/////////////////////////////////////////////////////
//Constructor/destructor
/////////////////////////////////////////////////////

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    //setup Qt UI
    setWindowIcon(QIcon ("://resources/icon.png"));
    ui->setupUi(this);

    //RJG - Add a few buttons to toolbar
    ui->mainToolBar->setIconSize(QSize(25,25));

    startButton = new QAction(QIcon(QPixmap(":/darkstyle/icon_play_button_green.png")), QString("Start"), this);
    startButton->setEnabled(true);
    ui->mainToolBar->addAction(startButton);ui->mainToolBar->addSeparator();
    QObject::connect(startButton, SIGNAL(triggered()), this, SLOT(on_actionStart_triggered()));

    stopButton = new QAction(QIcon(QPixmap(":/darkstyle/icon_stop_button_red.png")), QString("Stop"), this);
    stopButton->setEnabled(false);
    ui->mainToolBar->addAction(stopButton);ui->mainToolBar->addSeparator();
    QObject::connect(stopButton, SIGNAL(triggered()), this, SLOT(on_actionStop_triggered()));

    //RJG - add spacer
    QWidget* spacer = new QWidget();
    spacer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
    ui->mainToolBar->addWidget(spacer);

    aboutButton = new QAction(QIcon(QPixmap(":/darkstyle/icon_about_button.png")), QString("About"), this);
    ui->mainToolBar->addAction(aboutButton);
    stopButton->setEnabled(true);
    QObject::connect(aboutButton, SIGNAL (triggered()), this, SLOT (on_actionAbout_triggered()));

    //setup an action group to make tree export actions mutually exclusive
    alignmentGroup = new QActionGroup(this);
    alignmentGroup->addAction(ui->actionNewick);
    alignmentGroup->addAction(ui->actionPhyloXML);
    alignmentGroup->addAction(ui->actionTNT_tre);
    alignmentGroup->addAction(ui->actionTNT_MrBayes);
    ui->actionNewick->setChecked(true);

    //create a simulation opject
    TheSim=new Simulation;
    gotsomedata=false;

    setupGraphs();

    showMaximized();
}

MainWindow::~MainWindow()
{
    delete ui;
}

/////////////////////////////////////////////////////
//Methods associated with graphing system and tables
/////////////////////////////////////////////////////

void MainWindow::addGraph(int index, QString colour, int treemode)
{
    ui->widget->addGraph();
    ui->widget->graph(index)->setName(TheSim->modeToString(treemode));
    QPen gpen;
    gpen.setColor(QColor(colour));
    gpen.setWidth(4);
    ui->widget->graph(index)->setPen(gpen);
    graphindices[treemode]=index;
}

void MainWindow::setupGraphs()
{
    //setup QCustomPlot graph
    ui->widget->clearGraphs();
    graphindices.clear();
    for (int i=0; i<=TREE_MODE_MAX; i++) graphindices.append(-1); //empty index list
    int graphcount=0;
    if (ui->action_UseRDT->isChecked()) addGraph(graphcount++,"blue",TREE_MODE_RDT);
    if (ui->action_UseFDT->isChecked()) addGraph(graphcount++,"green",TREE_MODE_FDT);
    if (ui->action_UseIDT->isChecked()) addGraph(graphcount++,"red",TREE_MODE_IDT);
    QVector<qreal> dashes;
    qreal space = 8;
    dashes << 10 << space << 5 << space;

    // give the axes some labels:
    ui->widget->xAxis->setLabel("Species in genera");
    ui->widget->yAxis->setLabel("Frequency (%)");

    if (ui->actionLogged->isChecked())
    {
        ui->widget->xAxis->setScaleType(QCPAxis::stLogarithmic);
        ui->widget->xAxis->setScaleLogBase(10);
    }
    else
    {
        ui->widget->xAxis->setScaleType(QCPAxis::stLinear);
    }
    ui->widget->legend->setVisible(gotsomedata);
    ui->widget->setVisible(gotsomedata); //becomes visible when we get some data
    if (gotsomedata) plotcounts(&TheSim->counts,false);
}

void MainWindow::outputproportionalTable(QString name,QList<int> *data,QList<int> *satdata)
{
    QString s;
    QTextStream out(&s);
    std::sort(satdata->begin(),satdata->end()); //sort tree sizes
    out<<"\n\n------------------------------\nGenera sizes from "<<name<<" as binned percentages of total species in tree\n";
    out<<"------------------------------\nBetween,Count\n";
    csvdatabins[0].append(name+" proportional genus sizes");
    int base=0;
    for (int i=0; i<=PROPORTIONAL_BINS; i++)
    {
        double top=(double)(base+(100/PROPORTIONAL_BINS))-0.00001;
        if (base<100)
        {
                out<<base<<"%-"<<top<<"%,"<<data->at(i)<<"\n";
        }
        else
        {
            out<<"100%,"<<data->at(i)<<"\n";
        }
        csvdatabins[i+1].append(QString("%1").arg(data->at(i)));
        base+=(100/PROPORTIONAL_BINS);
    }
    out<<"----------------------------"<<"\n";
    out<<"Sizes of saturated (100%) trees:\n";
    for (int i=0; i<satdata->count(); i++)
    out<<satdata->at(i)<<"\n";
    out<<"----------------------------"<<"\n";
    QString treesizelist;
    for (int i=0; i<(satdata->count()-1); i++)
    treesizelist+=QString("%1~").arg(satdata->at(i));
    if (satdata->count()!=0) treesizelist+=QString("%1").arg(satdata->last());

    csvdatabins.last().append(treesizelist);
    logtext(s);

}

void MainWindow::outputTable(QString name, QVector<double> x, QVector<double>y, QHash<int,int>*count)
{
    QString s;
    QTextStream out(&s);
    out<<"\n\n------------------------------\nGenera sizes from "<<name<<"\n";
    bool nodata=false;
    if (count==0)
        nodata=true;
    else
       if (count->count()==0) nodata=true;

    //csv headers
    csvdata[0].append(name+"_freq");
    csvdata[0].append(name+"_perc");

    if (nodata)
        out<<"**** NO DATA ****\n";
    else
    {
        out<<"------------------------------\nGeneraSize,Frequency,Percentage\n";
        for (int i=0; i<x.length(); i++)
        {
            out<<x[i]<<","<<count->value(i+1,0)<<","<<y[i]<<"%\n";
            csvdata[i+1].append(QString("%1").arg(count->value(i+1,0)));
            csvdata[i+1].append(QString("%1%").arg(y[i]));
        }
    }
    out<<"----------------------------"<<"\n";
    logtext(s);

    //and - if csv data is turned out - put it in there too


}

void MainWindow::proportionaltables(QList<QList<int>*> *data,QList<QList<int>*> *satdata)
{
    //do table output
    if (ui->actionOutput_Genus_Size_to_Tree_Size_Tables->isChecked())
    {
        //set up csv headers for bins

        csvdatabins.clear();
        QStringList headers;
        headers.append("Bin >=");
        headers.append("Bin <");
        csvdatabins.append(headers);


        int base=0;
        for (int i=0; i<=PROPORTIONAL_BINS; i++)
        {
            double top=(double)(base+(100/PROPORTIONAL_BINS));
            if (base<100)
            {
                QStringList l;
                l.append(QString("%1%").arg(base));
                l.append(QString("%1%").arg(top));
                csvdatabins.append(l);
            }
            else
            {
                QStringList l;
                l.append(QString("100%"));
                l.append(QString("-"));
                csvdatabins.append(l);
            }
            base+=(100/PROPORTIONAL_BINS);
        }
        QStringList l;
        l.append("");
        l.append("100% tree sizes");
        csvdatabins.append(l);
        for (int i=1; i<=TREE_MODE_MAX; i++)
            if (graphindices[i]!=-1) outputproportionalTable(TheSim->modeToString(i),data->at(i),satdata->at(i));
    }

    if (getTableFileNeeded())
    {
        //write the csv
        QFile f(getCombinedTablesFileName());
        if (f.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            QTextStream out(&f);
            for (int i=0; i<csvdata.count(); i++)
            {
                out<<csvdata[i].join(",")<<"\n";
            }
            if (ui->actionOutput_Genus_Size_to_Tree_Size_Tables->isChecked())
            {
                out<<"\n-----genus size proportion tables-----\n";
                for (int i=0; i<csvdatabins.count(); i++)
                {
                    out<<csvdatabins[i].join(",")<<"\n";
                }
            }
            f.close();
        }
        else
            logtext(QString("Could not open combined tables file '%1' for writing").arg(f.fileName()));

    }

}

void MainWindow::plotcounts(QList<QHash<int,int>*> *data, bool showtable)
{
    //Updates the QCustomPlot chart.
    //Passed hash tables in a QList
        //positions are TREE_MODE codes
        //int contents is graph index for this code
        //For hash tables, key is number of genera, and whose value is count of occurrences of that number
    //showtable - table output now

    //expects the QList to have TREE_MODE_MAX + 1 entries - though they may be blank of course

    if (ui->actionOutputFrequencyTables->isChecked()==false) showtable=false;
    int maxcount=0;
    double maxy=0;

    //work out maximum count
    for (int ii=1; ii<data->length(); ii++)
        foreach(int i, data->at(ii)->keys())
            if (i>maxcount)maxcount=i;

    if (maxcount==0) return;

    //convert hash tables to lists of vectors, scaling properly and working out maxy
    QList<QVector<double> > datax, datay;
    for (int i=1; i<data->length(); i++) //skipping 0
    {
        QVector<double> x(maxcount),y(maxcount);

        double toty=0;
        QHash<int,int> *count=data->at(i);
        for (int ii=1; ii<=maxcount; ii++)
        {
            toty+=(y[ii-1]=count->value(ii,0));
            x[ii-1]=ii;
        }

        for (int ii=0; ii<maxcount; ii++)
        {
            y[ii]/=toty;
            y[ii]*=100;
            if (y[ii]>maxy) maxy=y[ii];
        }
        datax.append(x);
        datay.append(y);
    }

    csvdata.clear(); //empty any old csv data
    //set up correct number of lists, 0 is header

    QStringList l;
    l.append("Genus_size");
    csvdata.append(l);
    for (int i=1; i<=maxcount; i++)
    {
        QStringList l;
        l.append(QString("%1").arg(i));
        csvdata.append(l);
    }

    //do table output
    if (showtable)
    {
        for (int i=1; i<=TREE_MODE_MAX; i++)
            if (graphindices[i]!=-1) outputTable(TheSim->modeToString(i),datax[i-1],datay[i-1],data->value(i,(QHash<int,int>*)0));
    }

    //send data to graphs
    for (int i=1; i<=TREE_MODE_MAX; i++)
    {
        if (graphindices[i]!=-1) //-1 means graph not requested
            ui->widget->graph(graphindices[i])->setData(datax[i-1],datay[i-1]);
    }

    ui->widget->xAxis->setRange(1, maxcount);
    if (maxy>95) maxy=95;
    ui->widget->yAxis->setRange(0, maxy+5);
    ui->widget->legend->setVisible(true);
    ui->widget->setVisible(true);
    ui->widget->replot();
    gotsomedata=true;
}

void MainWindow::outputmaxgenussizefile(QVector< QList<maxgenusdatapoint> > *data)
{
    if (ui->actionOutput_max_genus_sizes_per_tree->isChecked())
    {
        QString s=TheSim->filepath+"/"+ui->stub->text()+"_maxgenussizedata.csv";
        QFile f(s);
        if (f.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            QTextStream out(&f);
            int maxcount=0;
            out<<"Iteration,tree_size";
            for (int i=1; i<=TREE_MODE_MAX; i++)
            {
                if ((*data)[i].count()>0) // there is data
                {
                    out<<","<<TheSim->modeToShortString(i);
                    maxcount=(*data)[i].count(); //should all be the same I hope
                    qDebug()<<maxcount;
                }
            }
            out<<"\n";

            for (int j=0; j<maxcount; j++)
            {
                out<<j;
                bool first=true;
                for (int i=1; i<=TREE_MODE_MAX; i++)
                {
                    if ((*data)[i].count()>0) // there is data
                    {
                        if (first)
                        {
                            out<<","<<(*data)[i].at(j).treesize;
                            first=false;
                        }
                        out<<","<<(*data)[i].at(j).maxgenussize;
                    }
                }
                out<<"\n";
            }
            f.close();
        }
        else logtext(QString("Couldn't open file '%1' for writing").arg(s));

    }
}

/////////////////////////////////////////////////////
//Other public output methods
/////////////////////////////////////////////////////


void MainWindow::logtext(QString text)
{
    //add string to log window
    ui->textEdit->append(text);
}

void MainWindow::setProgress(int current, int max)
{
    //set progress bar
    ui->progressBar->setMinimum(0);
    ui->progressBar->setMaximum(max);
    ui->progressBar->setValue(current);
    qApp->processEvents();
}


/////////////////////////////////////////////////////
//File menu handlers
/////////////////////////////////////////////////////

void MainWindow::on_actionExit_triggered()
{
    TheSim->stop();
    qApp->quit();
    exit(0);
}

void MainWindow::on_actionStart_triggered()
{
    startButton->setEnabled(false);
    stopButton->setEnabled(true);
    TheSim->run(this);
}

void MainWindow::on_actionStop_triggered()
{
    startButton->setEnabled(true);
    stopButton->setEnabled(false);
    TheSim->stop();
}

void MainWindow::on_actionChart_to_PDF_triggered()
{
    QString filename=QFileDialog::getSaveFileName(this,"Output file name",TheSim->filepath,"PDF Files (*.pdf)");
    if (filename=="") return;
    ui->widget->savePdf(filename);
}



/////////////////////////////////////////////////////
//Tree export menu handlers
/////////////////////////////////////////////////////

void MainWindow::on_actionSet_Export_Folder_2_triggered()
{
   QString d = QFileDialog::getExistingDirectory(this,"Folder for output");
   if (d!="") TheSim->filepath=d;
}


/////////////////////////////////////////////////////
//Mode menu handlers
/////////////////////////////////////////////////////

void MainWindow::on_actionLogged_triggered()
{
    //change scale between log/linear
    if (ui->actionLogged->isChecked())
    {
        ui->widget->xAxis->setScaleType(QCPAxis::stLogarithmic);
        ui->widget->xAxis->setScaleLogBase(10);
        ui->widget->replot();
    }
    else
    {
        ui->widget->xAxis->setScaleType(QCPAxis::stLinear);
        ui->widget->replot();
    }
}



/////////////////////////////////////////////////////
//Tree file functions
/////////////////////////////////////////////////////


void MainWindow::log_identical_genomes(QString text)
//part of the specific leaf count code. Used to log a per-tree identical genome count to a file
{
    QString s=TheSim->filepath+"/"+ui->stub->text()+"_identicals.txt";
    QFile f(s);
    if (f.open(QIODevice::Append | QIODevice::Text))
    {
        f.write((QString("%1\n").arg(text)).toLatin1());
        f.close();
    }

}

void MainWindow::do_with_matrix_trees(int iter, Lineage *rootlineage)
{
    //public, called by simulation. Exports the with matrix trees to file - not done
    //per taxonomy, so in different function
    if (ui->actionExport_Trees->isChecked()==false || ui->actionTNT_MrBayes->isChecked()==false) return; //nothing to do


    QString s=TheSim->filepath+"/"+ui->stub->text()+"_batch.nex";
    QFile f(s);
    if (f.open(QIODevice::Append))
    {
        f.write(TheSim->dumpnex_with_matrix(rootlineage,iter).toLatin1());
        f.close();
    }
    s=TheSim->filepath+"/"+ui->stub->text()+"_batch.tnt";
    QFile f2(s);
    if (f2.open(QIODevice::Append))
    {
        f2.write(TheSim->dumptnt_with_matrix(rootlineage,iter).toLatin1());
        f2.close();
    }
    return;

}


void MainWindow::do_trees(int mode, int iter, Lineage *rootlineage)
{
    //public, called by simulation. Export trees to file

    if (getTaxonomyTypeInUse(mode)==false || ui->actionExport_Trees->isChecked()==false) return; //nothing to do

    QString extension=".tre";
    QString iters;
    iters.sprintf("%06d",iter);
    if (ui->actionTNT_MrBayes->isChecked()) iters.sprintf("%06d",iter-1);

    if (ui->actionNewick->isChecked()) extension=".nwk";
    if (ui->actionPhyloXML->isChecked()) extension=".phyloxml";
    if (ui->actionTNT_MrBayes->isChecked()) extension=".nex";

    QString s=TheSim->filepath+"/"+ui->stub->text()+"_"+TheSim->modeToShortString(mode)+"_"+iters+extension;
    QFile f(s);
    if (f.open(QIODevice::WriteOnly))
    {
        if (ui->actionNewick->isChecked()) f.write(TheSim->dumpnewick(rootlineage).toLatin1());
        if (ui->actionPhyloXML->isChecked()) f.write(TheSim->dumpphyloxml(rootlineage).toLatin1());
        if (ui->actionTNT_tre->isChecked()) f.write(TheSim->dumptnt(rootlineage).toLatin1());
        if (ui->actionTNT_MrBayes->isChecked()) f.write(TheSim->dump_nex_alone(rootlineage).toLatin1());
        f.close();
    }
}


/////////////////////////////////////////////////////
//Gets for UI parameters / settings
/////////////////////////////////////////////////////

double MainWindow::getchanceextinct()
{
   return ui->extinctchance->value();
}

double MainWindow::getchancespeciate()
{
   return ui->speciatechance->value();
}


int MainWindow::getabsthreshold()
{
   return ui->threshold_abs->value();
}


double MainWindow::getRDTthreshold()
{
   return ui->RDT_threshold->value();
}

int MainWindow::getIDTthreshold()
{
   return ui->threshold_idt->value();
}

int MainWindow::getgenerations()
{
   return ui->generations->value();
}
int MainWindow::getiterations()
{
   return ui->time_increments->value();
}
int MainWindow::getmaxleafcount()
{
   return ui->maxleaves->value();
}

double MainWindow::getchancemutate()
{
   return ui->mutatechance->value();
}

int MainWindow::get_precise_leaf_count()
{
    return ui->precise_leaf_count->value();
}

bool MainWindow::getTableFileNeeded()
{
    return ui->actionCombine_Tables_into_File->isChecked();
}

QString MainWindow::getCombinedTablesFileName()
{

    if (TheSim->filepath.length()==0) return "";
    else
    {
        QString fname=ui->stub->text();
        fname+="_tables.csv";
        return TheSim->filepath+"/"+fname;
    }
}

bool MainWindow::getTaxonomyTypeInUse(int mode)
{
   if (mode==TREE_MODE_FDT) return ui->action_UseFDT->isChecked();
   if (mode==TREE_MODE_RDT) return ui->action_UseRDT->isChecked();
   if (mode==TREE_MODE_IDT) return ui->action_UseIDT->isChecked();
   if (mode==TREE_MODE_TNTMB) return true; //bit of a bodge - but this always writes trees
   return false;
}

bool MainWindow::throw_away_on_leaf_limit()
{
    if (ui->actionThrow_away_trees_hitting_leaf_limit->isChecked()) return true; else return false;
}

bool MainWindow::per_tree_output()
{
    if (ui->actionOutput_result_per_tree->isChecked()) return true; else return false;
}

bool MainWindow::correct_number_trees_wanted()
{
    if (ui->actionContinue_until_correct_number_of_trees->isChecked()) return true; else return false;
}


void MainWindow::on_action_UseFDT_triggered()
{
    setupGraphs();
}

void MainWindow::on_action_UseIDT_triggered()
{
    setupGraphs();
}

void MainWindow::on_action_UseRDT_triggered()
{
    setupGraphs();
}

void MainWindow::on_actionAbout_triggered()
{
    QMessageBox::about(this,"About MBL2017 v2","MBL2017 v2 (with additions for character evolution). coding Mark Sutton, m.sutton@ic.ac.uk<br />Concepts and algorithms Mark Sutton, Julia Sigwart.");
}
