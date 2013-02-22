Disicc::Application.routes.draw do
  devise_for :users, :path => "user", :path_names => { :sign_in => 'login', :sign_out => 'logout', :password => 'secret', :confirmation => 'verification', :unlock => 'unblock', :registration => 'register', :sign_up => 'sign_up' }
  resources :users
  resources :page
  resources :sequences do
    member do
      get :run_disorder
      get :disorder_consensus
      get :download_disorder
    end
  end
  resources :a_asequencs
  
  resources :alignments do
    member do
      get :display_annotated_alignment
      get :display_compensatory_annotated_alignment
      get :display_xdet_annotated_alignment
      get :display_caps_annotated_alignment
      get :display_disorder_annotated_alignment
      get :calculate_intraresidue_consensus
      get :calculate_intraresidue_consensus_threaded
      get :calculate_disorder_consensus_threaded
      get :percent_identities
      get :compensatory_brief_report
      get :disorder_brief_report
    end
    collection do
      get :upload
      post :pre_process_fasta_file 
      post :process_fasta_file 
      post :complete_process_fasta_file 
      post :process_fasta_file_and_save_sequences
    end
     #match 'alignments/:id', :to => 'catalog#display_annotated_alignment'
  end
 
  resources :alignment_positions
  resources :caps
  resources :classifications
  resources :conseqs
  resources :disorder_type
  resources :disorders
  resources :disorder_values
  resources :intra_residue_contacts
  resources :xdets
  # The priority is based upon order of creation:
  # first created -> highest priority.

  # Sample of regular route:
  #   match 'products/:id' => 'catalog#view'
  # Keep in mind you can assign values other than :controller and :action

  # Sample of named route:
  #   match 'products/:id/purchase' => 'catalog#purchase', :as => :purchase
  # This route can be invoked with purchase_url(:id => product.id)

  # Sample resource route (maps HTTP verbs to controller actions automatically):
  #   resources :products

  # Sample resource route with options:
  #   resources :products do
  #     member do
  #       get 'short'
  #       post 'toggle'
  #     end
  #
  #     collection do
  #       get 'sold'
  #     end
  #   end

  # Sample resource route with sub-resources:
  #   resources :products do
  #     resources :comments, :sales
  #     resource :seller
  #   end

  # Sample resource route with more complex sub-resources
  #   resources :products do
  #     resources :comments
  #     resources :sales do
  #       get 'recent', :on => :collection
  #     end
  #   end

  # Sample resource route within a namespace:
  #   namespace :admin do
  #     # Directs /admin/products/* to Admin::ProductsController
  #     # (app/controllers/admin/products_controller.rb)
  #     resources :products
  #   end

  # You can have the root of your site routed with "root"
  # just remember to delete public/index.html.
  root :to => "page#index"

  # See how all your routes lay out with "rake routes"

  # This is a legacy wild controller route that's not recommended for RESTful applications.
  # Note: This route will make all actions in every controller accessible via GET requests.
  # match ':controller(/:action(/:id(.:format)))'
end
