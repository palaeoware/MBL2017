#include "simulation.h"
#include "lineage.h"
#include "mainwindow.h"
#include <QString>
#include <QDebug>
#include <QTime>
#include <QFile>
#include <QFileInfo>
#include <QMutableHashIterator>
#include <QMutableLinkedListIterator>
#include <QInputDialog>
#include "genus.h"
#include "math.h"
#include <QMessageBox>

/////////////////////////////////////////////////////
//Simulation class - single-instance class to manage simulation
/////////////////////////////////////////////////////


/////////////////////////////////////////////////////
//Global data declarations
/////////////////////////////////////////////////////

Simulation *TheSimGlobal;
QHash<qint64,Genus*> genera;
int CHANCE_EXTINCT;
int CHANCE_SPECIATE;
double CHANCE_SPECIATE_DOUBLE, CHANCE_EXTINCT_DOUBLE;

int ABS_THRESHOLD;
double RDT_THRESHOLD;
int IDT_THRESHOLD;
quint32 CHANCE_MUTATE;
double dCHANCE_MUTATE;

quint64 rtot;
quint64 rcount;
qint64 genusnumberidt;

quint32 tweakers[32]; //used for modifying characters during evolution
int bitcounts[65536]; //bitcount array - for each possible 16 bit word, how many bits are on?

Lineage *dummy_parameter_lineage;


/////////////////////////////////////////////////////
//Constructor/Destructor and run
/////////////////////////////////////////////////////

Simulation::Simulation()
{
    rootspecies=0;
    crownroot=0;
    TheSimGlobal=this;
    filepath="c:/";  //CHANGE THIS FOR NON-WINDOWS SYSTEMS!

    //set up list of data arrays - one for each possible type and a dummy for 0 (unclassified)
    for (int i=0; i<=TREE_MODE_MAX; i++) counts.append(new QHash<int,int>);
    for (int i=0; i<=TREE_MODE_MAX; i++) proportional_counts.append(new QList<int>);
    for (int i=0; i<=TREE_MODE_MAX; i++) saturated_tree_sizes.append(new QList<int>);

}

Simulation::~Simulation()
{

}

void Simulation::run(MainWindow *mainwin)
{

    //performs a simulation run

    //keep pointer to main window instance
    mw=mainwin;

    stopflag=false;
    rfilepoint=0;

    //setup random stuff
    rtot=0;
    rcount=0;
    read65536numbers();

    rpoint=0;

    //do tweakers
    tweakers[0]=1;
    for (int n=1; n<32; n++) tweakers[n]=tweakers[n-1]*2;

    //and bitcounts array
    for (qint32 n=0; n<65536; n++)
    {
            qint32 count=0;
            for (int m=0; m<16; m++) if ((n & tweakers[m])!=0) ++count;  // count the bits
            bitcounts[n]=count;
    }

    //set up max genus size per tree
    QVector< QList<maxgenusdatapoint> > maxgenussize; //List of lists, done per taxonomy type. MS 32 bits - tree size. LS 32 bits - max genus size

    for (int i=0; i<=TREE_MODE_MAX; i++)
    {
        QList<maxgenusdatapoint> l;
        maxgenussize.append(l);
    }

    //clear any data from previous runs
    if (genera.count()>0)
        qDeleteAll(genera);
    genera.clear();
    for (int i=0; i<counts.length(); i++) counts.at(i)->clear(); //clear count data
    for (int i=0; i<proportional_counts.length(); i++)
    {
        proportional_counts.at(i)->clear();
        for (int j=0; j<=PROPORTIONAL_BINS; j++) //yes, <= is correct - one extra bin for 100%
            proportional_counts.at(i)->append(0); //set up empty bins
    }

    for (int i=0; i<saturated_tree_sizes.length(); i++) saturated_tree_sizes.at(i)->clear();
   //get settings from UI, in int and double types - int type is in range 0-65535
    CHANCE_EXTINCT=(int)(mw->getchanceextinct()*65535.0);
    CHANCE_SPECIATE=(int)(mw->getchancespeciate()*65535.0);
    CHANCE_EXTINCT_DOUBLE=mw->getchanceextinct();
    CHANCE_SPECIATE_DOUBLE=mw->getchancespeciate();

    ABS_THRESHOLD=mw->getabsthreshold();
    RDT_THRESHOLD=mw->getRDTthreshold();
    IDT_THRESHOLD=mw->getIDTthreshold();

    //Get mutation chance
    dCHANCE_MUTATE=mw->getchancemutate();
    CHANCE_MUTATE=(quint32)(dCHANCE_MUTATE*(double)65536.0*(double)65536.0);

    int generations=mw->getgenerations(); //this is number of trees to run
    int iterations=mw->getiterations(); //time iterations per tree
    int leafcountmax=mw->getmaxleafcount();

    int preciseleafcount=mw->get_precise_leaf_count();

    //run counts/totals
    int totalca=0; //count alive
    int totalce=0; //count extinct
    int totalcb=0; //count branched (speciated)
    int actualtreecount=0;
    int actualiterations=0;

    if (preciseleafcount) mw->log_identical_genomes("tree,identical_genomes");

    //main loop - iterate for specified number of trees
    for (int i=0; i<generations; i++)
    {
        actualiterations++;
        QString so;
        QTextStream sout(&so);
        sout<<"Iteration "<<i<<" actualiterations "<<actualiterations<<" actualtreecount "<<actualtreecount;
        mw->logtext(so+"\n");

        dummy_parameter_lineage=0; //important - constructor will try to check it!
        quint32 dummycharacters[CHARACTER_WORDS];
        randomcharacters(dummycharacters);
        dummy_parameter_lineage=new Lineage(dummycharacters,(Lineage *)0,0,0); //run with single species that lacks a parent, created at time step 0

        if (mw->correct_number_trees_wanted()) i=0; //if iterating to correct number of trees ensure loop does not terminate on for statement

        //progress bar
        if (mw->correct_number_trees_wanted())
            mw->setProgress(actualtreecount,generations);
        else
            mw->setProgress(i,generations);

        //run data for run
        leafcount=1;
        nextgenusnumber=1;
        nextid=0;
        currenttime=0;
        if (rootspecies!=0) delete rootspecies; //delete any existing root species from previous iteration of main loop
        if (crownroot!=0) delete crownroot; //delete any existing fossil free tree from previous iteration of main loop
        crownroot=0;

        //create a root species
        quint32 startcharacters[CHARACTER_WORDS];
        randomcharacters(startcharacters);
        rootspecies=new Lineage(startcharacters,(Lineage *)0,0,0); //run with single species that lacks a parent, created at time step 0


        //Run tree for correct number of iterations
        int j;
        for (j=0; j<iterations; j++)
        {
            increment();
            if (leafcount>leafcountmax) break;
        }

        if (dummy_parameter_lineage) delete dummy_parameter_lineage;


        //check - did it terminate because it hit leaf limit?
        if (j!=iterations)
        {
            //yes - log this
            QString so;
            QTextStream sout(&so);
            sout<<"Leaf limit hit: Tree "<<i<<" ran for "<<j<<" iterations, hit "<<leafcount<<" leaves.";
            mw->logtext(so+"\n");
            if (mw->throw_away_on_leaf_limit())
            {
                if (mw->per_tree_output()) //if tree is to be discarded AND per tree output is on, log the discard
                {
                    QString out;
                    out.sprintf("%d: Discarded (too many leaves)",actualiterations);
                    mw->logtext(out);
                }
                continue; //skip to next iteration of main loop
            }
        }

        if (preciseleafcount!=0) //cull to precise size if possible
        {
            crownroot=rootspecies->strip_extinct((Lineage *)0); //removes extinct tips
            if (crownroot!=0)
            {
                crownroot->cull_dead_branches(); //merge single branch nodes

                Lineage *preciseleafcountroot = crownroot->find_clade_with_precise_size(preciseleafcount);
                if (preciseleafcountroot)
                {
                     mw->logtext("Found a precise clade");
                    //found one!

                    //check - how many uninformative characters - discard if this is not 0
                    int unin=0;
                    for (int i=0; i<CHARACTER_WORDS; i++)
                    {
                        quint32 mask=1;
                        for (int j=0; j<32; j++)
                        {
                            int z=preciseleafcountroot->count_zeros(i,mask);
                            int o=preciseleafcountroot->count_ones(i,mask);
                            //qDebug()<<"0"<<z<<"1"<<o;
                            if (z>=(preciseleafcount-1)) unin++; //if all zeroes except 1 or 0 1's - uninformative
                            if (o>=(preciseleafcount-1)) unin++; //if all ones except 0 or 1 0's - uninformative
                            mask*=2;
                        }
                    }

                    if (unin>0)
                        mainwin->logtext("Discarding tree - contains one or more uninformative characters");
                    else
                    {
                        actualtreecount++;

                        mw->logtext("OK, found one");
                        //find out if there are any identical genomes in there
                        QSet<QString> genomes;
                        preciseleafcountroot->addstringgenomestoset(&genomes);
                        QString identical_genome_string=QString("%1,%2").arg(actualtreecount-1).arg(2*(preciseleafcount-genomes.count()));
                        mw->log_identical_genomes(identical_genome_string);

                        nextsimpleIDnumber=0;
                        preciseleafcountroot->dosimplenumbers(); //do simple numbering scheme
                        mw->do_with_matrix_trees(actualtreecount, preciseleafcountroot);
                        mw->do_trees(TREE_MODE_TNTMB,actualtreecount,preciseleafcountroot);
                        mw->logtext(QString("Found subtree with length %1\n").arg(preciseleafcountroot->count_alive()));
                    }
                }
                else  mw->logtext("COULD NOT findd a precise clade");

            }
        }
        else
        {
            //normal mode (no precise count specified - as old code)
            //get count of how many exant, extinct and branched lineages descend from the root species
            int ca=rootspecies->count_alive();
            int ce=rootspecies->count_extinct();
            int cb=rootspecies->count_branched();
            totalca+=ca; //add to totals
            totalce+=ce;
            totalcb+=cb;

            //log leaf count of tree if per-tree output is on
            if (mw->per_tree_output())
            {
                QString out;
                out.sprintf("Iteration %d: extant leaves at end of run: %d",actualiterations,ca);
                mw->logtext(out);
            }

            if (ca>0)
            {
                //If there is a tree, i.e. we have any extant descendents of the root species

                actualtreecount++;
                //Produce version with extinct stripped
                crownroot=rootspecies->strip_extinct((Lineage *)0); //removes extinct tips
                if (crownroot!=0) crownroot->cull_dead_branches(); //merge single branch nodes

                nextsimpleIDnumber=0;
                crownroot->dosimplenumbers(); //do simple numbering scheme
                mw->do_with_matrix_trees(i, crownroot);

                //Now go through the various taxonomy modes

                //1. TREE_MODE_UNCLASSIFIED
                mw->do_trees(TREE_MODE_UNCLASSIFIED,i,rootspecies); //handle the unclassified tree writing

                leafcount=ca;

                //2. TREE_MODE_RDT
                if (mw->getTaxonomyTypeInUse(TREE_MODE_RDT))
                {
                    //clear genus list
                    if (genera.count()>0)  qDeleteAll(genera); //delete all existing genera data structures
                    genera.clear();
                    //no need to clear classifications here - crownroot not yet used
                    RDT_genera(); //do RDT classification

                    int max=genera_data_report(TREE_MODE_RDT);
                    maxgenusdatapoint dp;
                    dp.maxgenussize=max;
                    dp.treesize=leafcount;
                    maxgenussize[TREE_MODE_RDT].append(dp);

                    mw->do_trees(TREE_MODE_RDT,i,crownroot);
                }

                //3. TREE_MODE_FDT
                if (mw->getTaxonomyTypeInUse(TREE_MODE_FDT))
                {
                    //clear genus list
                    if (genera.count()>0)  qDeleteAll(genera); //delete all existing genera data structures
                    genera.clear();
                    crownroot->resetGenusLabels();
                    crownroot->genuslabels_distance(0,j); //do FDT classification
                    int max=genera_data_report(TREE_MODE_FDT);
                    maxgenusdatapoint dp;
                    dp.maxgenussize=max;
                    dp.treesize=leafcount;
                    maxgenussize[TREE_MODE_FDT].append(dp);
                    mw->do_trees(TREE_MODE_FDT,i,crownroot);
                }

                //4. TREE_MODE_IDT
                if (mw->getTaxonomyTypeInUse(TREE_MODE_IDT))
                {
                    //clear genus list
                    if (genera.count()>0)  qDeleteAll(genera); //delete all existing genera data structures
                    genera.clear();
                    genusnumberidt=1;
                    crownroot->resetGenusLabels();
                    crownroot->genuslabels_IDT(genusnumberidt,j); //do IDT classification
                    int max=genera_data_report(TREE_MODE_IDT);
                    maxgenusdatapoint dp;
                    dp.maxgenussize=max;
                    dp.treesize=leafcount;
                    maxgenussize[TREE_MODE_IDT].append(dp);
                    mw->do_trees(TREE_MODE_IDT,i,crownroot);
                }

                mw->plotcounts(&counts,false); //plot the results
            }
        }

        //should we stop?
        if (mw->correct_number_trees_wanted() && actualtreecount==generations) stopflag=true;
        if (stopflag) {generations=i; break;}

    }

    //set progress bar to complete
    mw->plotcounts(&counts,true); //plot the results with tables too - resets CSV
    mw->proportionaltables(&proportional_counts,&saturated_tree_sizes);  // - also outputs csv
    mw->setProgress(generations,generations);

    mw->outputmaxgenussizefile(&maxgenussize);

}


/////////////////////////////////////////////////////
//Random number methods
/////////////////////////////////////////////////////

quint32 Simulation::GetRandom16()
{
    qint32 r=(quint32)((GetRandom() / 65536));
    return r;
}

quint32 Simulation::GetRandom()
{
    if (rpoint<65535)
    return randoms[rpoint++];
    else
    {
        read65536numbers();
        rpoint=0;
        return randoms[rpoint++];
    }
}

void Simulation::read65536numbers()
{
    //read more random numbers from file
    int c[10];
    rfilepoint+=65536*4;
#ifdef Q_OS_MAC  //should be outside package
    QFile f("../../../randomnumbers.dat");
#else //win32 or linux
    QFile f("randomnumbers.dat");
#endif
    QFileInfo f2(f.fileName());

    //qDebug()<<f2.absoluteFilePath();
    if (f.open(QIODevice::ReadOnly)==false)
    {
        QMessageBox::warning(mw,"No random numbers!","No random number file found - see readme file. This should be named 'randomnumbers.dat', and placed on path (in "+f2.absoluteFilePath()+"). MBL will now exit - please sort out this file and restart.");
        exit(0);
    }
    if ((f.size()-rfilepoint)<65536*4) rfilepoint=qrand();
    f.seek(rfilepoint);
    f.read((char *)(&(randoms[0])),65536*4);
    f.close();

    for (int i=0; i<10; i++) c[i]=0;

    qint64 tot=0;
    for (int i=0; i<65536; i++)
    {
        double r=(double)randoms[i];
        r/=65536;
        r/=65536;
        r*=10;
        int ri=(int)r;
        c[ri]++;
        tot+=(qint64)randoms[i];
    }
}



/////////////////////////////////////////////////////
//Daughter methods for run
/////////////////////////////////////////////////////

void Simulation::increment()
{
    //one iteration of rootspecies
    currenttime++;
    if (rootspecies!=0)
    {
        rootspecies->iterate(currenttime); //will chain-iterate them all
    }
}

void Simulation::RDT_genera()
{
    if (crownroot==0) return;
    //calculate genera using recursive RDT rules

    if (genera.count()>0) qDeleteAll(genera);
    genera.clear();

    QList<Lineage *>extantlist;
    crownroot->getextantlist(&extantlist);

    qint64 gnumber=1;

    for (int i=0; i<extantlist.count(); i++)
    {
        Lineage *thisone=extantlist[i];
        if (thisone->parent_lineage) //exclude root!
        {
            Lineage *sister=thisone->getsister();

            if (sister==0) {continue;}

            if (sister->time_split==-1 && sister->genusnumber==0) //not split and must still be alive - so simple sister - though exclude if already labelled
            {
                //simple monophyletic pair. Create as a genus.
                Genus *g= new Genus;
                g->id=gnumber;
                g->species.append(thisone);
                g->species.append(sister);
                thisone->genusnumber=gnumber;
                sister->genusnumber=gnumber;
                g->rootnode=thisone->parent_lineage;
                gnumber++;
                genera.insert(g->id,g);

                bool continueloop=true;
                do
                {
                    //look at sister clade to this genus
                    sister=g->rootnode->getsister();


                    if (sister==0)
                    {
                        //this happens if we try to expand out past root node. We can't - just move on
                        break;
                    }

                    int agegenusMRCA=currenttime-g->rootnode->time_split;
                    int agesistergenusMRCA=currenttime-g->rootnode->time_created;
                    if ((int)((double)agegenusMRCA/RDT_THRESHOLD)>=agesistergenusMRCA)
                    {
                        //passed the threshold cut off rule for depth - might be incorporatable
                        if (sister->RDT_check(currenttime))
                        {
                            sister->RDT_incorporate(g);
                            g->rootnode=g->rootnode->parent_lineage;
                            if (g->rootnode==crownroot) continueloop=false;
                        }
                        else
                            continueloop=false;
                    }
                    else
                        continueloop=false;
                }
                while (continueloop);
            }
        }
    }

    //deal with remaining singletons
    for (int i=0; i<extantlist.count(); i++)
    {
        if (extantlist[i]->genusnumber==0)
        {
            Genus *g= new Genus;
            g->id=gnumber;
            g->species.append(extantlist[i]);
            extantlist[i]->genusnumber=gnumber;
            gnumber++;
            genera.insert(g->id,g);
        }
    }


}


int Simulation::get_mrca_age(QList<Lineage *> *terminals)
{
    //takes a list of nodes, and returns the age of the mrca.
    if (terminals->length()<2) qDebug("ERROR in get_mrca_age");
    Lineage *mrca=get_mrca(terminals->at(0),terminals->at(1));
    for (int i=2; i<terminals->length(); i++)
    {
        if (mrca->isThisADescendent(terminals->at(i))) continue; //skip
        else
            mrca=get_mrca(mrca,terminals->at(i));
    }
    return currenttime-mrca->time_split;
}

Lineage * Simulation::get_mrca(Lineage *l0, Lineage *l1)
{
    //determine most recent common ancestor (MRCA) for two lineages
    QList <Lineage *> l0ancestors, l1ancestors;

    while(l0->parent_lineage)
    {
        l0ancestors.append(l0->parent_lineage);
        l0=l0->parent_lineage;
    };
    while(l1->parent_lineage)
    {
        l1ancestors.append(l1->parent_lineage);
        l1=l1->parent_lineage;
    };

    //got lists, find first in L0 list that is also in L1 list

    foreach(Lineage *l, l0ancestors)
    {
        if (l1ancestors.contains(l)) return l;
    }

    //shouldn't get here!
    qDebug()<<"INTERNAL ERROR: No MRCA found in Simulation::get_mrca";
    return (Lineage *)0;
}

int Simulation::genera_data_report(int mode)
//returns maximum genus size found
{

    int max=0;
    //add stuff to counts hash table for graphing and tables
    QHash<int,int> *thiscounts=counts[mode];
    foreach (Genus *g, genera)
    {
        int oldval=thiscounts->value(g->species.count(),0);
        (*thiscounts)[g->species.count()]=oldval+1;
    }

    QList<int> *thispcounts=proportional_counts[mode];
    foreach (Genus *g, genera)
    {
        if (g->species.count()>max) max=g->species.count();
        double proportion = ((double)g->species.count())/(double)leafcount;
        int bin=(int) (proportion * (100/PROPORTIONAL_BINS));
        if (g->species.count()==leafcount)
        {
            saturated_tree_sizes.at(mode)->append(leafcount);
            bin=PROPORTIONAL_BINS; //make sure 100% is 100%
        }
        (*thispcounts)[bin]=thispcounts->at(bin)+1;
    }

    return max;
}

int Simulation::distancebetween(quint32 chars1[], quint32 chars2[])
//works out distance between a pair of genomes
{
    int total=0;
    for (int i=0; i<CHARACTER_WORDS; i++)
    {
        quint32 c1=chars1[i];
        quint32 c2=chars2[i];
        quint32 diffs = c1 ^ c2; //XOR the two to compare
        total+= bitcounts[diffs/(quint32)65536] +  bitcounts[diffs & (quint32)65535];
    }
    return total;
}

/////////////////////////////////////////////////////
//Miscellaneous
/////////////////////////////////////////////////////

void Simulation::stop()
{
    stopflag=true; //triggers a stop at end of current tree
}


QString Simulation::modeToString(int mode)
{
    if (mode==TREE_MODE_FDT) return "Fixed Depth Taxonomy (FDT)";
    if (mode==TREE_MODE_RDT) return "Relative Top-Down Taxonomy (RDT)";
    if (mode==TREE_MODE_UNCLASSIFIED) return "No Taxonomy";
    if (mode==TREE_MODE_IDT) return "Internal Distance Taxonomy (IDT)";
    if (mode==TREE_MODE_TNTMB) return "TNT/MB .nex output";

    return "error-typenotfound";
}

QString Simulation::modeToShortString(int mode)
{
    if (mode==TREE_MODE_FDT) return "FDT";
    if (mode==TREE_MODE_RDT) return "RDT";
    if (mode==TREE_MODE_UNCLASSIFIED) return "NoTax";
    if (mode==TREE_MODE_IDT) return "IDT";
    if (mode==TREE_MODE_TNTMB) return "TNTMB";

    return "error-typenotfound";
}


QString Simulation::dumptnt_with_matrix(Lineage *rootlineage, int iter)
{
    QString ret;

    //work out proper numbering system

    QString n = QString("%1").arg(nextsimpleIDnumber-1);
    int numberofdigits=n.length();
    QString zerostring;
    zerostring.sprintf(QString("%0"+QString("%1").arg(numberofdigits)+"d").toUtf8(),0);


    ret+="xread\n";
    ret+="'Written on "+QDateTime::currentDateTime().toString("ddd MMM d hh:mm:ss yyyy")+"'\n";
    ret+=QString("%1 %2\n").arg(32*CHARACTER_WORDS).arg(nextsimpleIDnumber);
    ret+=rootlineage->getcharactermatrix(false);
    ret+=";\n\n";
    ret+="hold 10000\n";
    ret+="mult=tbr replic 1000 hold 1000;\n";
    ret+=QString("export-treesim_%1_POUT.tre;\nxwipe;\n\n").arg(iter-1);
    return ret;
}

QString Simulation::dumpnex_with_matrix(Lineage *rootlineage, int iter)
{
    QString ret;
    QString n = QString("%1").arg(nextsimpleIDnumber);
    int numberofdigits=n.length();
    QString zerostring;
    zerostring.sprintf(QString("%0"+QString("%1").arg(numberofdigits)+"d").toUtf8(),0);

    ret+="#NEXUS\n";
    ret+="[Written on "+QDateTime::currentDateTime().toString("ddd MMM d hh:mm:ss yyyy")+"]\n";
    ret+=QString("\nBegin data;\nDimensions ntax=%2 nchar=%1;\n").arg(32*CHARACTER_WORDS).arg(nextsimpleIDnumber);
    ret+="Format datatype=standard missing=? gap=-;\n\nMatrix\n";
    ret+=rootlineage->getcharactermatrix(false); //not root - i.e. no outgroup
    ret+=";\nEND;\n\n";
    ret+="begin mrbayes;\n";
    ret+="        set autoclose=yes nowarn=yes;\n";
    ret+="        lset coding=variable rates=gamma;\n";
    ret+="        \n";
    ret+="        [mcmc settings]\n";
    ret+="        set usebeagle=yes beagledevice=CPU beagleprecision=double beaglescaling=always beaglesse=no beagleopenmp=no;\n";
    ret+="        \n";
    ret+="        mcmcp temp=0.1 nchain=4 samplefreq=200 printfr=1000 nruns=2 append=no;\n";
    ret+="        \n";
    ret+=QString("        mcmcp filename=%1;\n").arg(iter-1);
    ret+="        \n";
    ret+="        mcmc ngen=2000000;\n";
    ret+="        \n";
    ret+="        sumt relburnin =yes burninfrac = 0.25; [set relative burnin to 25% for consensus tree]\n";
    ret+="        \n";
    ret+="        sump relburnin =yes burninfrac = 0.25; [set relative burnin to 25% for tree probabilities]\n";
    ret+="        \n";
    ret+="        end;\n";
    ret+="        \n";
    return ret;
}

QString Simulation::dumpnewick(Lineage *rootlineage)
{
    return rootlineage->newickstring();
}

QString Simulation::dumptnt(Lineage *rootlineage)
{
    return "tread 'tree dumped from TreeModel in TNT format'\n"+rootlineage->tntstring()+"\nproc-;";
}

QString Simulation::dump_nex_alone(Lineage *rootlineage)
{
    QString ret;
    ret="#NEXUS\n";
    ret+="Begin trees;\n";
    ret+="    Translate\n";
    for (int i=0; i<nextsimpleIDnumber; i++)
    {
        QString comma;
        if (i!=(nextsimpleIDnumber-1)) comma=",";
        ret+=QString("%1 Species_%2%3\n").arg(i+1).arg(i).arg(comma);
    }
    ret+=";\n\n\n";

    QString tree=rootlineage->numbertree(1);
    //need to strip off last :number

    //qDebug()<<"BEFORE"<<tree;
    for (int i=tree.length()-1; i>0; i--)
    {
        if (tree[i]==QChar(':'))
        {
            tree=tree.left(i);
            break;
        }
    }
    //qDebug()<<"AFTER"<<tree;
    ret+="tree tree1 =[&U]"+tree+";\n";
    ret+="END;\n";
    return ret;
}

QString Simulation::dumpphyloxml(Lineage *rootlineage)
{
    return "<?xml version=\"1.0\" encoding=\"UTF-8\"?><phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd\" xmlns=\"http://www.phyloxml.org\"> <phylogeny rooted=\"true\" rerootable=\"true\">"
            + rootlineage->xmlstring()
            + "</phylogeny></phyloxml>\n";
}

void Simulation::randomcharacters(quint32 *chars)
//generates randomized characters for initial lineage
{

    for (int i=0; i<CHARACTER_WORDS; i++)
    {
        chars[i]=GetRandom();
    }
}
