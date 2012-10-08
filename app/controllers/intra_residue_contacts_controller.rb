class IntraResidueContactsController < ApplicationController
  # GET /intra_residue_contacts
  # GET /intra_residue_contacts.xml
  def index
    @intra_residue_contacts = IntraResidueContact.all

    respond_to do |format|
      format.html # index.html.erb
      format.xml  { render :xml => @intra_residue_contacts }
    end
  end

  # GET /intra_residue_contacts/1
  # GET /intra_residue_contacts/1.xml
  def show
    @intra_residue_contact = IntraResidueContact.find(params[:id])

    respond_to do |format|
      format.html # show.html.erb
      format.xml  { render :xml => @intra_residue_contact }
    end
  end

  # GET /intra_residue_contacts/new
  # GET /intra_residue_contacts/new.xml
  def new
    @intra_residue_contact = IntraResidueContact.new

    respond_to do |format|
      format.html # new.html.erb
      format.xml  { render :xml => @intra_residue_contact }
    end
  end

  # GET /intra_residue_contacts/1/edit
  def edit
    @intra_residue_contact = IntraResidueContact.find(params[:id])
  end

  # POST /intra_residue_contacts
  # POST /intra_residue_contacts.xml
  def create
    @intra_residue_contact = IntraResidueContact.new(params[:intra_residue_contact])

    respond_to do |format|
      if @intra_residue_contact.save
        format.html { redirect_to(@intra_residue_contact, :notice => 'Intra residue contact was successfully created.') }
        format.xml  { render :xml => @intra_residue_contact, :status => :created, :location => @intra_residue_contact }
      else
        format.html { render :action => "new" }
        format.xml  { render :xml => @intra_residue_contact.errors, :status => :unprocessable_entity }
      end
    end
  end

  # PUT /intra_residue_contacts/1
  # PUT /intra_residue_contacts/1.xml
  def update
    @intra_residue_contact = IntraResidueContact.find(params[:id])

    respond_to do |format|
      if @intra_residue_contact.update_attributes(params[:intra_residue_contact])
        format.html { redirect_to(@intra_residue_contact, :notice => 'Intra residue contact was successfully updated.') }
        format.xml  { head :ok }
      else
        format.html { render :action => "edit" }
        format.xml  { render :xml => @intra_residue_contact.errors, :status => :unprocessable_entity }
      end
    end
  end

  # DELETE /intra_residue_contacts/1
  # DELETE /intra_residue_contacts/1.xml
  def destroy
    @intra_residue_contact = IntraResidueContact.find(params[:id])
    @intra_residue_contact.destroy

    respond_to do |format|
      format.html { redirect_to(intra_residue_contacts_url) }
      format.xml  { head :ok }
    end
  end
end
